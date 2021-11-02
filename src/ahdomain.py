# Configuation class for ah domains in HOOMD
# XXX: Maybe make a dataclass in the future
import hoomd
import hoomd.md as md
import gsd.hoomd

import itertools
import os
import random
import sys
import yaml

import numpy as np
from common import *

class ahdomain(object):
    def __init__(self, bead_size, yaml_file):
        # AH-domain parameters
        self.nah              = np.int32(np.float64(yaml_file['ah_domain']['n_ah']))
        self.polymer_type     = yaml_file['ah_domain']['polymer_type']
        self.gamma            = np.float64(yaml_file['ah_domain']['gamma'])
        self.nbeads           = np.int32(np.float64(yaml_file['ah_domain']['nbeads']))
        self.nrepeat          = np.int32(np.float64(yaml_file['ah_domain']['nrepeat']))
        self.mass             = np.float64(yaml_file['ah_domain']['mass'])
        self.A                = np.float64(yaml_file['ah_domain']['A'])
        self.Bsurface         = np.float64(yaml_file['ah_domain']['B_surface'])
        self.Bintermediate    = np.float64(yaml_file['ah_domain']['B_intermediate'])
        self.Bdeep            = np.float64(yaml_file['ah_domain']['B_deep'])
        self.Bself            = np.float64(yaml_file['ah_domain']['B_self'])
        self.kbond            = np.float64(yaml_file['ah_domain']['kbond'])
        self.kbend            = np.float64(yaml_file['ah_domain']['kbend'])
        self.is_init          = yaml_file['ah_domain']['is_init']

        # Set constants we will need
        self.bead_size = 0.75*bead_size
        self.mass_per_bead = self.mass / self.nbeads
        self.r0 = self.bead_size
        self.rc = 2.0*self.r0
        self.rbond = self.r0

    # Initialize the AH domains
    def InitAH(self, snap):
        if self.is_init:
            print(f"Already think AH is initialized")
        else:
            self.CreateAH(snap)
        self.is_init = True

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
        # Make sure we avoid problems with the box overlap and the domains
        rboxxy = snap.configuration.box[0]/2.0 - self.bead_size*self.nbeads
        rboxz = snap.configuration.box[2]/2.0 - self.bead_size*self.nbeads
        # Loop over number of AH domains
        for ahdx in range(self.nah):
            # Put in some overlap logic
            while True:
                overlap = 0
                # Assign the location of the AH filament
                xrand = random.uniform(-1.0*rboxxy, rboxxy)
                yrand = random.uniform(-1.0*rboxxy, rboxxy)
                zrand = random.uniform(self.bead_size*self.nbeads, rboxz)
                # Assign a random unit vector
                v = generate_random_unit_vector(3)

                ah_start_x = np.array([xrand, yrand, zrand])
                ndx = 0
                ah_code = False
                for idx in range(ah_start_nx + self.nbeads*ahdx, ah_start_nx + self.nbeads*(ahdx+1)):
                    if ndx % self.nrepeat == 0:
                        ah_code = not ah_code

                    # Calculate the position
                    pos = ah_start_x + v*self.bead_size*ndx

                    # Check for overlap with other AH domains
                    for jdx in range(ah_start_nx, ah_start_nx + self.nbeads*ahdx):
                        r1 = snap.particles.position[jdx]
                        overlap += np.int32(sphere_overlap(self.bead_size, pos, r1))
                        if overlap != 0:
                            break

                    if overlap != 0:
                        break

                    snap.particles.position[idx] = pos
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

                # Bail for the while loop
                if overlap == 0:
                    break
       
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

