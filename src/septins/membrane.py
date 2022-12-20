# XXX: Copyright info here

"""Membrane class for septin project"""

import hoomd
import hoomd.md as md
import gsd.hoomd

import itertools
import os
import random
import sys
import yaml

import numpy as np

class Membrane(object):
    def __init__(self, verbose, bead_size, yaml_file):
        self.verbose = verbose
        if self.verbose: print(f"Membrame::__init__")
        # Lipid parameters
        self.nbeads    = np.int32(np.float64(yaml_file['membrane']['nbeads']))
        self.mdopc     = np.float64(yaml_file['membrane']['mass'])
        self.gamma     = np.float64(yaml_file['membrane']['gamma'])
        self.A         = np.float64(yaml_file['membrane']['A'])
        self.B         = np.float64(yaml_file['membrane']['B'])
        self.kbond     = np.float64(yaml_file['membrane']['kbond'])
        self.kbend     = np.float64(yaml_file['membrane']['kbend'])
        self.is_init   = yaml_file['membrane']['is_init']
        self.ngrid     = np.int32(np.float64(yaml_file['membrane']['ngrid']))
        self.lbox      = np.float64(yaml_file['simulation']['lbox'])
        self.ndirect   = np.int32(np.float64(yaml_file['membrane']['ndirect']))
        self.bead_size = bead_size

        # Head beads are slightly smaller than the normal beads
        self.r0 = self.bead_size
        self.rc = 2.0 * self.r0

        # Also account for a slightly different diffusion for them
        if 'gamma_head' in yaml_file['membrane']:
            self.gamma_head = np.float64(yaml_file['membrane']['gamma_head'])
        else:
            self.gamma_head = self.gamma * 0.75
   
        # Special choices for the head bead size
        self.rhh = 0.75 * self.r0
        self.Ahh = 0.75 * self.A
        
        # Bond distance and length, set to the distance between molecules
        self.rbond = self.r0
        
        # Lipid mass in g/mol 
        self.lipidmass = self.mdopc
        self.lipidmass_per_bead = self.lipidmass / self.nbeads

        # Initial configuration information for analysis
        self.analysis_init = False
        self.getTypebyName = {}
        self.getNamebyType = {}
        self.nmembrane = self.nbeads * 2 * self.ngrid * self.ngrid
        self.nlipids = 2 * self.ngrid * self.ngrid
        self.nlipids_per_leaflet = self.nlipids / 2

        if self.verbose: print(f"Membrame::__init__ return")

    # Init any information that required a snapshot
    def InitMembrane(self, snap):
        if self.verbose: print(f"Membrane::InitMembrane")
        if self.is_init:
            print(f"Creating membrane from file")
            for itype,ntype in enumerate(snap.particles.types):
                self.getTypebyName[ntype] = itype
        else:
            print(f"Creating membrane from scratch!")
            self.getTypebyName = {'H': 0, 'I': 1, 'T': 2}
            self.CreateMembrane(snap)
        self.is_init = True
        if self.verbose: print(f"Membrane::InitMembrane return")

    # Create the membrane for our system
    def CreateMembrane(self, snap):
        if self.verbose: print(f"Membrane::CreateMembrane")
        # What is the number of replications we need in the box?
        box_extent = self.lbox
        linear_extent = box_extent/self.ngrid
        fudge_size = 0.2 # Prevent particles on the boundaries because reasons
        snap.configuration.box = [linear_extent, linear_extent, 2.0*((self.bead_size/2.0) + (self.nbeads-1)*self.bead_size)+fudge_size, 0, 0, 0]

        # Create the number of bead types we need for just membranes
        snap.particles.N = 2*self.nbeads
        snap.particles.types = ['H', 'I', 'T']

        # The number of bonds we have is determined by nbeads, as is the number of angle potentials
        snap.bonds.N = 2*(self.nbeads - 1)
        snap.angles.N = 2*(self.nbeads - 2)
        snap.bonds.types = ['lipidbond']
        snap.bonds.typeid[:] = 0
        snap.angles.types = ['lipidbend']
        snap.angles.typeid[:] = 0

        # We have two leaflets, assign each with a zleaf loop
        # Start idx at 0, this will increment for every change!
        idx = 0
        for zidx in [1, -1]:
            for i in range(self.nbeads):
                snap.particles.position[idx] = [0.0, 0.0, zidx*( (self.bead_size/2.0) + (self.nbeads-1-i)*self.bead_size )]
        
                # Set the typeid in an ugly if else block, but easiest way to control the double Tail portion
                if i == 0:
                    snap.particles.typeid[idx] = self.getTypebyName['H']
                elif i == 1:
                    snap.particles.typeid[idx] = self.getTypebyName['I']
                else:
                    snap.particles.typeid[idx] = self.getTypebyName['T']
        
                # Set the mass for the beads to be the same
                snap.particles.mass[idx] = self.lipidmass_per_bead
        
                idx += 1
        
        # Sadly, hardcoding the bond types for each number of beads seems best
        if self.nbeads == 3:
            snap.bonds.group[:] = [[0, 1], [1, 2], [3, 4], [4, 5]]
            snap.angles.group[:] = [[0, 1, 2], [3, 4, 5]]
        elif self.nbeads == 4:
            snap.bonds.group[:] = [[0, 1], [1, 2], [2, 3], [4, 5], [5, 6], [6, 7]]
            snap.angles.group[:] = [[0, 1, 2], [1, 2, 3], [4, 5, 6], [5, 6, 7]]
        else:
            print("Only used for 3,4-bead models, exiting", file=sys.stderr)
        
        # Replicate the snapshot on a lattice
        snap.replicate(self.ngrid, self.ngrid, 1)
        
        # Create the box size of the system to be something reasonable
        snap.configuration.box = [box_extent, box_extent, box_extent, 0, 0, 0]

        # XXX: FIXME: Dump a lammps configuration at the same time!
        self.CreateLAMMPS(snap)
        if self.verbose: print(f"Membrane::CreateMembrane return")

    def CreateLAMMPS(self, snap):
        if self.verbose: print(f"Membrane::CreateLAMMPS")

        # Get the  number of atoms
        natoms = snap.particles.N
        ntypes = len(snap.particles.types)
        nbonds = snap.bonds.N
        nbondstypes = self.nbeads - 1
        nangles = snap.angles.N
        nanglestypes = self.nbeads - 2
        lx = snap.configuration.box[0]
        ly = snap.configuration.box[1]
        lz = snap.configuration.box[2]

        lammps_str = (
            f'Autogenerated from dragonfruit\n'
            f'\n'
            f'{natoms} atoms\n'
            f'{ntypes} atom types\n'
            f'{nbonds} bonds\n'
            f'{nbondstypes} bond types\n'
            f'{nangles} angles\n'
            f'{nanglestypes} angle types\n'
            f'\n'
            f'  {-lx/2} +{lx/2} xlo xhi\n'
            f'  {-ly/2} +{ly/2} ylo yhi\n'
            f'  {-lz/2} +{lz/2} zlo zhi\n'
            f'\n'
            f'Masses\n\n'
        )

        # Loop through types to get masses
        for itype in range(ntypes):
            lammps_str += f'  {itype+1} {self.lipidmass_per_bead}\n'
        lammps_str += f'\n'

        # Loop over all atoms in system, their molecule, and their type
        lammps_str += f'Atoms\n\n'
        imolecule = -1 
        itype = 0
        for iatom in range(natoms):
            ibead = iatom % self.nbeads
            if ibead == 0: imolecule += 1
            # The logic for the tail type is just ugly, brute force
            if ibead == 0:
                itype = 0
            elif ibead == 1:
                itype = 1
            else:
                itype = 2

            lammps_str += f'  {iatom+1}  {imolecule+1}  {itype+1} {snap.particles.position[iatom][0]:4.2f} {snap.particles.position[iatom][1]:4.2f} {snap.particles.position[iatom][2]:4.2f}\n'

        lammps_str += f'\n'

        # Loop over the bonds in the system to write them too
        lammps_str += f'Bonds\n\n'
        for ibond in range(nbonds):
            bondtype = ibond % nbondstypes
            bondgroup = snap.bonds.group[ibond]
            lammps_str += f'  {ibond+1}  {bondtype+1}  {bondgroup[0]+1}  {bondgroup[1]+1}\n'

        lammps_str += f'\n'
        
        # Loop over the angles and generate
        lammps_str += f'Angles\n\n'
        for iangle in range(nangles):
            angletype = iangle % nanglestypes
            anglegroup = snap.angles.group[iangle]
            lammps_str += f'  {iangle+1}  {angletype+1}  {anglegroup[0]+1}  {anglegroup[1]+1}  {anglegroup[2]+1}\n'

        lammps_str += f'\n'

        lammps_filename = f'{self.ngrid}x{self.ngrid}_membrane.lammps_config'
        with open(lammps_filename, 'w') as stream:
            stream.write(lammps_str)

        if self.verbose: print(f"Membrane::CreateLAMMPS return")

    def PrintInformation(self, snap):
        # Print out the relevant information about the lipids
        print(f"--------")
        print(f"Lipid information")
        print(f"N beads in membrane     = {self.nmembrane}")
        print(f"N lipids                = {self.nlipids}")
        print(f"Lipid mass (amu)        = {self.lipidmass}")
        print(f"  Lipid mass per bead   = {self.lipidmass_per_bead}")
        print(f"Gamma (m/tau)           = {self.gamma}")
        print(f"  Gamma head (m/tau)    = {self.gamma_head}")
        print(f"A (kBT)                 = {self.A}")
        print(f"B (kBT)                 = {self.B}")
        print(f"R (sigma)               = {self.r0}")
        print(f"RC (sigma)              = {self.rc}")
        print(f"(Ahh (kBT))             = {self.Ahh}")
        print(f"(Rhh (sigma))           = {self.rhh}")
        print(f"Bond r (sigma)          = {self.rbond}")
        print(f"Bond k (kBT/sigma^2)    = {self.kbond}")
        print(f"Bend k (kBT/rad^2)      = {self.kbend}")

        lbox = snap.configuration.box[0]
        self.area_per_lipid = (lbox * lbox) / (self.nlipids_per_leaflet)
        print(f"--------")
        print(f"Derived membrane values")
        print(f"  Number of lipid patches:              {self.nlipids}")
        print(f"  Number of lipid patches per leaflat:  {self.nlipids_per_leaflet}")
        print(f"  Area per lipid:                       {self.area_per_lipid}")

        print(f"--------")
        print(f"Analysis parameters")
        print(f"  Ndirect:                              {self.ndirect}")

    def ConfigureAnalysis(self, snap):
        r""" Configure the analysis for the trajectory
        """
        if self.verbose: print(f"Membrane::ConfigureAnalysis")

        if self.analysis_init:
            print(f"Membrane: analysis already set up!")
            return

        # Get the TypeID by name or by type
        for idx,val in enumerate(snap.particles.types):
            self.getTypebyName[val] = idx
            self.getNamebyType[idx] = val

        if self.verbose: print(f"Membrane::ConfigureAnalysis return")

    def GetLeafletIndices(self, snap, ptype):
        r""" Get the indices for the lipids, upper and lower leaflets
        """
        idx = np.where(snap.particles.typeid == self.getTypebyName[ptype])[0]
        leaf1_idx = idx[:len(idx)//2]
        leaf2_idx = idx[len(idx)//2:]

        return [idx, leaf1_idx, leaf2_idx]

