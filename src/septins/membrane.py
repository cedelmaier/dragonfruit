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
    def __init__(self, bead_size, yaml_file):
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
        self.bead_size = bead_size

        # Head beads are slightly smaller than the normal beads
        self.r0 = self.bead_size
        self.rc = 2.0 * self.r0
   
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
        self.getTypeIDbyName = {}
        self.getNamebyTypeID = {}
        self.nmembrane = self.nbeads * 2 * self.ngrid * self.ngrid
        self.nlipids = 2 * self.ngrid * self.ngrid

    # Init any information that required a snapshot
    def InitMembrane(self, snap):
        if self.is_init:
            print(f"Creating membrane from file")
            self.getTypeByName = {}
            for itype,ntype in enumerate(snap.particles.types):
                self.getTypeByName[ntype] = itype
        else:
            print(f"Creating membrane from scratch!")
            self.getTypeByName = {'H': 0, 'I': 1, 'T': 2}
            self.CreateMembrane(snap)
        self.is_init = True

    # Create the membrane for our system
    def CreateMembrane(self, snap):
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
                    snap.particles.typeid[idx] = self.getTypeByName['H']
                elif i == 1:
                    snap.particles.typeid[idx] = self.getTypeByName['I']
                else:
                    snap.particles.typeid[idx] = self.getTypeByName['T']
        
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

    def PrintInformation(self, snap):
        # Print out the relevant information about the lipids
        print(f"--------")
        print(f"Lipid information")
        print(f"N beads in membrane     = {self.nmembrane}")
        print(f"N lipids                = {self.nlipids}")
        print(f"Lipid mass (amu)        = {self.lipidmass}")
        print(f"  Lipid mass per bead   = {self.lipidmass_per_bead}")
        print(f"Gamma (m/tau)           = {self.gamma}")
        print(f"A (kBT)                 = {self.A}")
        print(f"B (kBT)                 = {self.B}")
        print(f"R (sigma)               = {self.r0}")
        print(f"RC (sigma)              = {self.rc}")
        print(f"(Ahh (kBT))             = {self.Ahh}")
        print(f"(Rhh (sigma))           = {self.rhh}")
        print(f"Bond r (sigma)          = {self.rbond}")
        print(f"Bond k (kBT/sigma^2)    = {self.kbond}")
        print(f"Bend k (kBT/rad^2)      = {self.kbend}")

        self.nlipids_per_leafleat = self.nlipids / 2
        lbox = snap.configuration.box[0]
        self.area_per_lipid = (lbox * lbox) / (self.nlipids_per_leafleat)
        print(f"--------")
        print(f"Derived membrane values")
        print(f"  Number of lipid patches:              {self.nlipids}")
        print(f"  Number of lipid patches per leafleat: {self.nlipids_per_leafleat}")
        print(f"  Area per lipid:                       {self.area_per_lipid}")

    def ConfigureAnalysis(self, snap):
        r""" Configure the analysis for the trajectory
        """
        print(f"Membrane::ConfigureAnalysis")

        if self.analysis_init:
            print(f"Membrane: analysis already set up!")
            return

        # Get the TypeID by name or by type
        for idx,val in enumerate(snap.particles.types):
            self.getTypeIDbyName[val] = idx
            self.getNamebyTypeID[idx] = val

        print(f"Membrane::ConfigureAnalysis return")

    def GetLeafletIndices(self, snap, ptype):
        r""" Get the indices for the lipids, upper and lower leaflets
        """
        idx = np.where(snap.particles.typeid == self.getTypeIDbyName[ptype])[0]
        leaf1_idx = idx[:len(idx)//2]
        leaf2_idx = idx[len(idx)//2:]

        return [idx, leaf1_idx, leaf2_idx]

