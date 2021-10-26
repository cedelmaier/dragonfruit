# Configuration class for Lipid systems (extensible in future)
import hoomd
import hoomd.md as md
import gsd.hoomd

import itertools
import os
import random
import sys
import yaml

import numpy as np

# The configurator really exists as simulation to convert between HOOMD
# and our other simulatino variables

# Simulation parameters for a membrane object
class Membrane(object):
    def __init__(self, bead_size, yaml_file):
        # Lipid parameters
        self.nbeads    = np.int32(np.float32(yaml_file['membrane']['nbeads']))
        self.mdopc     = np.float32(yaml_file['membrane']['mass'])
        self.gamma     = np.float32(yaml_file['membrane']['gamma'])
        self.A         = np.float32(yaml_file['membrane']['A'])
        self.B         = np.float32(yaml_file['membrane']['B'])
        self.kbond     = np.float32(yaml_file['membrane']['kbond'])
        self.kbend     = np.float32(yaml_file['membrane']['kbend'])
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

    # Create the membrane for our system
    def CreateMembrane(self, snap, ngrid, lbox):
        # What is the number of replications we need in the box?
        box_extent = lbox
        linear_extent = box_extent/ngrid
        fudge_size = 0.2 # Prevent particles on the boundaries because reasons

        # Create the number of bead types we need for just membranes
        snap.particles.N = 2*self.nbeads
        snap.particles.types = ['H', 'I', 'T']
        self.getTypeByName = {'H': 0, 'I': 1, 'T': 2}

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
        snap.replicate(ngrid, ngrid, 1)
        
        # Create the box size of the system to be something reasonable
        snap.configuration.box = [box_extent, box_extent, box_extent, 0, 0, 0]

        return True
       
    def PrintInformation(self, snap):
        # Print out the relevant information about the lipids
        print(f"--------")
        print(f"Lipid information")
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

        # Write out any interesting derived information
        self.nlipids = 0
        for idx in range(snap.particles.N):
            if snap.particles.typeid[idx] == 0:
                self.nlipids += 1
        self.nlipids_per_leafleat = self.nlipids / 2
        lbox = snap.configuration.box[0]
        self.area_per_lipid = (lbox * lbox) / (self.nlipids_per_leafleat)
        print(f"--------")
        print(f"Derived membrane values")
        print(f"  Number of lipid patches:              {self.nlipids}")
        print(f"  Number of lipid patches per leafleat: {self.nlipids_per_leafleat}")
        print(f"  Area per lipid:                       {self.area_per_lipid}")

# Information on a block copolymer model
class AHBlockCopolymer(object):
    def __init__(self, bead_size, yaml_file):
        # AH-domain parameters
        self.nah              = np.int32(np.float32(yaml_file['ah_domain']['n_ah']))
        self.polymer_type     = yaml_file['ah_domain']['polymer_type']
        self.gamma            = np.float32(yaml_file['ah_domain']['gamma'])
        self.nbeads           = np.int32(np.float32(yaml_file['ah_domain']['nbeads']))
        self.nrepeat          = np.int32(np.float32(yaml_file['ah_domain']['nrepeat']))
        self.mass             = np.float32(yaml_file['ah_domain']['mass'])
        self.A                = np.float32(yaml_file['ah_domain']['A'])
        self.Bsurface         = np.float32(yaml_file['ah_domain']['B_surface'])
        self.Bintermediate    = np.float32(yaml_file['ah_domain']['B_intermediate'])
        self.Bdeep            = np.float32(yaml_file['ah_domain']['B_deep'])
        self.Bself            = np.float32(yaml_file['ah_domain']['B_self'])
        self.kbond            = np.float32(yaml_file['ah_domain']['kbond'])
        self.kbend            = np.float32(yaml_file['ah_domain']['kbend'])

        # Set constants we will need
        self.bead_size = 0.75*bead_size
        self.mass_per_bead = self.mass / self.nbeads
        self.r0 = self.bead_size
        self.rc = 2.0*self.r0
        self.rbond = self.r0

        # Set up the RNG
        random.seed(10)

    # Create some number of block copolyerms in the simulation
    def CreateAH(self, snap):
        # Early bail if no AH domains
        if self.nah == 0:
            return

        # We have additional types that we must setup before creating particles
        p_ntypes = len(snap.particles.types)
        snap.particles.types = snap.particles.types + ['AH1', 'AH2']
        self.getTypeByName = {'AH1': p_ntypes+0, 'AH2': p_ntypes+1}
        b_ntypes = len(snap.bonds.types)
        snap.bonds.types = snap.bonds.types + ['ahbond']
        self.getTypeByName['ahbond'] = b_ntypes
        a_ntypes = len(snap.angles.types)
        snap.angles.types = snap.angles.types + ['ahbend']
        self.getTypeByName['ahbend'] = a_ntypes

        # This is for the copolymer model, some linear number of bonds
        ah_start_nx = snap.particles.N
        snap.particles.N = snap.particles.N + (self.nbeads * self.nah)
        # Loop over number of AH domains
        for ahdx in range(self.nah):
            # Assign the location of the AH filament
            xrand = random.uniform(-1.0*snap.configuration.box[0]/2, 1.0*snap.configuration.box[0]/2)
            yrand = random.uniform(-1.0*snap.configuration.box[1]/2, 1.0*snap.configuration.box[1]/2)

            print(f"Inserting AH-domain start at ({xrand, yrand})")

            ah_start_x = np.array([xrand, yrand, snap.configuration.box[0]/4.0])
            ndx = 0
            ah_code = False
            for idx in range(ah_start_nx + self.nbeads*ahdx, ah_start_nx + self.nbeads*(ahdx+1)):
                if ndx % self.nrepeat == 0:
                    ah_code = not ah_code
                snap.particles.position[idx] = ah_start_x + np.array([0, (self.bead_size/2.0)*(idx - ah_start_nx), 0])
                if ah_code:
                    snap.particles.typeid[idx] = self.getTypeByName['AH1']
                else:
                    snap.particles.typeid[idx] = self.getTypeByName['AH2']
                snap.particles.mass[idx] = self.mass_per_bead

                ndx += 1

            # Now set up the bond information
            n_newbonds = self.nbeads - 1
            ah_start_bdx = snap.bonds.N
            bdx = ah_start_bdx
            snap.bonds.N = snap.bonds.N + n_newbonds
            for idx in range(ah_start_nx + self.nbeads*ahdx, ah_start_nx + self.nbeads*ahdx + self.nbeads-1):
                snap.bonds.typeid[bdx] = self.getTypeByName['ahbond']
                snap.bonds.group[bdx] = [idx, idx+1]
                bdx += 1

            # Now set up the angle information
            n_newangles = self.nbeads - 2
            ah_start_adx = snap.angles.N
            adx = ah_start_adx
            snap.angles.N = snap.angles.N + n_newangles
            for idx in range(ah_start_nx + self.nbeads*ahdx, ah_start_nx + self.nbeads*ahdx + self.nbeads - 2):
                snap.angles.typeid[adx] = self.getTypeByName['ahbend']
                snap.angles.group[adx] = [idx, idx+1, idx+2]
                adx += 1

       
    def PrintInformation(self, snap):
        # Print out AH information
        print(f"--------")
        print(f"AH information (copolymer model)")
        print(f"AH number                           = {self.nah}")
        print(f"AH bead size (sigma)                = {self.bead_size}")
        print(f"AH mass (amu)                       = {self.mass}")
        print(f"  AH mass per bead                  = {self.mass_per_bead}")
        print(f"Gamma (m/tau)                       = {self.gamma}")
        print(f"AH A (kBT / sigma)                  = {self.A}")
        print(f"AH B terms (kBT / sigma)")
        print(f"  B surface                         = {self.Bsurface}")
        print(f"  B intermediate                    = {self.Bintermediate}")
        print(f"  B deep                            = {self.Bdeep}")
        print(f"  B self                            = {self.Bself}")
        print(f"AH R (sigma)                        = {self.r0}")
        print(f"AH RC (sigma)                       = {self.rc}")
        print(f"AH r (sigma)                        = {self.rbond}")
        print(f"AH kbond (kBT/sigma^2)              = {self.kbond}")
        print(f"AH kbend (kBT/rad^2)                = {self.kbend}")
    

# Class definition
class Configurator(object):
    def __init__(self, opts):
        self.default_yaml = self.GetYamlDict(opts.default_file)

        self.SetDefaults()

    def SetDefaults(self):
        # Simulation parameters
        self.kT                 = np.float32(self.default_yaml['simulation']['kT'])
        self.bead_size          = np.float32(self.default_yaml['simulation']['bead_size'])
        self.deltatau           = np.float32(self.default_yaml['simulation']['deltatau'])
        self.nsteps             = np.int32(np.float32(self.default_yaml['simulation']['nsteps']))
        self.nwrite             = np.int32(np.float32(self.default_yaml['simulation']['nwrite']))
        self.ngrid              = np.int32(np.float32(self.default_yaml['simulation']['ngrid']))
        self.lbox               = np.float32(self.default_yaml['simulation']['lbox'])
        self.nseed              = np.int32(np.float32(self.default_yaml['simulation']['seed']))
        self.compute_mode       = self.default_yaml['simulation']['mode']
        self.trajectory_file    = self.default_yaml['simulation']['trajectory_file']
        self.init_type          = self.default_yaml['simulation']['init_type']
        self.is_membrane_init   = self.default_yaml['simulation']['is_membrane_init']
        self.is_ahdomain_init   = self.default_yaml['simulation']['is_ahdomain_init']
        self.integrator         = self.default_yaml['simulation']['integrator']

        # Check integrator settings to make sure it's okay!
        if (self.integrator != 'langevin' and self.integrator != 'NPT' and self.integrator != 'NPH'):
            print(f"Integrator {self.integrator} not used for this system, exiting!")
            sys.exit(1)

        # Unfortunately, now have some gymnastics for checking if something exists and setting
        if 'pdamp' in self.default_yaml['simulation']:
            self.pdamp = np.float32(self.default_yaml['simulation']['pdamp'])
        else:
            self.pdamp = None

        if 'tau' in self.default_yaml['simulation']:
            self.tau = np.float32(self.default_yaml['simulation']['tau'])
        else:
            self.tau = None

        if 'tauS' in self.default_yaml['simulation']:
            self.tauS = np.float32(self.default_yaml['simulation']['tauS'])
        else:
            self.tauS = None

        if self.init_type == 'read_gsd':
            self.init_filename = self.default_yaml['simulation']['init_filename']

        # Lipid parameters
        if 'membrane' in self.default_yaml:
            self.lipids = Membrane(self.bead_size, self.default_yaml)

        # AH parameters
        if 'ah_domain' in self.default_yaml:
            if self.default_yaml['ah_domain']['polymer_type'] == 'block_copolymer':
                self.ahdomain = AHBlockCopolymer(self.bead_size, self.default_yaml)
            else:
                print("Only block copolymer available for now, exiting!")
                sys.exit(1)
        else:
            print("No AH domain included, exiting!")
            sys.exit(1)

    # Print information
    def PrintInformation(self, snap):
        print(f"--------")
        print(f"System information")
        print(f"Compute mode            = {self.compute_mode}")
        print(f"Integrator              = {self.integrator}")
        print(f"Delta tau               = {self.deltatau}")
        print(f"Simulation time (tau)   = {self.deltatau*self.nsteps}")
        print(f"kBT                     = {self.kT}")
        print(f"seed                    = {self.nseed}")
        print(f"box                     = ({snap.configuration.box[0]}, {snap.configuration.box[1]}, {snap.configuration.box[2]})")
        # Configurable parameters
        if self.tau: print(f"tau                     = {self.tau}")
        if self.pdamp: print(f"pdamp                   = {self.pdamp}")
        if self.tauS: print(f"tauS                    = {self.tauS}")

        self.lipids.PrintInformation(snap)
        self.ahdomain.PrintInformation(snap)

    # Passthrough for creation options
    def CreateMembrane(self, snap):
        if not self.is_membrane_init:
            self.is_membrane_init = self.lipids.CreateMembrane(snap, self.ngrid, self.lbox)
    def CreateAH(self, snap):
        if not self.is_ahdomain_init:
            self.is_ahdomain_init = self.ahdomain.CreateAH(snap)

    # Get a YAML dictionary for the file
    def GetYamlDict(self, filename):
        file_dict = ''
        with open(filename, 'r') as stream: file_dict = yaml.safe_load(stream)
        return file_dict

