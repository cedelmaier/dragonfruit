# Configuration class for helix AH domains in HOOMD
import hoomd
import hoomd.md as md
import gsd.hoomd

import itertools
import os
import random
import sys
import yaml

import numpy as np
from scipy.optimize import fsolve
from ahdomain import ahdomain
from common import *

def helix_position(a, b, t):
    return np.array([a*np.cos(t), a*np.sin(t), b*t], dtype=np.float64)

def helix_space_func(t0, a, b, d):
    return 4.0 * np.square(a) * np.square(np.sin(t0/2.0)) + np.square(b*t0) - np.square(d)

class helix_ahdomain(ahdomain):
    def __init__(self, bead_size, yaml_file):
        super().__init__(bead_size, yaml_file)

        # Override the number of beads, etc, that we are using
        self.pitch                          = np.float64(yaml_file['ah_domain']['pitch'])
        self.radius                         = np.float64(yaml_file['ah_domain']['radius'])
        self.bead_size                      = np.float64(yaml_file['ah_domain']['bead_size'])
        self.translation_per_residue        = np.float64(yaml_file['ah_domain']['translation_per_residue'])
        self.naa                            = np.float64(yaml_file['ah_domain']['naa'])
        self.initial_arrangement            = yaml_file['ah_domain']['initial_arrangement']
        self.sequence                       = yaml_file['ah_domain']['sequence']
        self.initial_location               = yaml_file['ah_domain']['initial_location']

        # Derived quantities
        self.total_length = self.naa * self.translation_per_residue
        a = self.radius
        b = self.pitch / (2*np.pi)
        self.curvature = np.abs(self.radius)/(self.radius**2 + b**2)
        tspace_guess = 2.0 / np.sqrt(self.radius**2 + b**2) / self.curvature * np.arcsin(self.curvature*self.bead_size / 2.0)
        tspace_root = fsolve(helix_space_func, tspace_guess, args=(a, b, self.bead_size))
        self.tspace = tspace_root[0]

        self.nbeads = np.int64(self.total_length / (self.tspace * b))
        self.mass_per_bead = self.mass / self.nbeads
        self.r0 = self.bead_size
        self.rc = 2.0*self.r0
        self.rbond = self.r0
        self.beads_per_residue = self.nbeads / self.naa

        r01 = helix_position(self.radius, b, self.tspace) - helix_position(self.radius, b, 0.0)
        r12 = helix_position(self.radius, b, 2.0*self.tspace) - helix_position(self.radius, b, self.tspace)
        theta0 = np.arccos(np.dot(r01, r12)/np.linalg.norm(r01)/np.linalg.norm(r12))
        self.thetabend = np.pi - theta0

        # Make some calculations based on the hydrophobic nature of the sequence
        self.hydrophobic_dict = {'A': True,
                                 'I': True,
                                 'L': True,
                                 'M': True,
                                 'F': True,
                                 'V': True,
                                 'P': True,
                                 'G': True}

        if self.initial_location != "random" and self.initial_location != "upper_half":
            print(f"Initial location {self.initial_location} not recognized, exiting!")
            sys.exit(1)

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

        # This is for the helix model, some linear number of bonds
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
                if self.initial_location == "upper_half":
                    zrand = random.uniform(self.bead_size*self.nbeads, rboxz)
                elif self.initial_location == "random":
                    zrand = random.uniform(-1.0*rboxz, rboxz)
                # Assign a random unit vector
                v = generate_random_unit_vector(3)

                ah_start_x = np.array([xrand, yrand, zrand])
                ndx = 0 # Place in current AH domain
                residue_idx = 0 # Current residue number
                #ah_code = False
                for idx in range(ah_start_nx + self.nbeads*ahdx, ah_start_nx + self.nbeads*(ahdx+1)):
                    #if ndx % self.nrepeat == 0:
                    #    ah_code = not ah_code

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

                    # Figure out if the bead is AH1 or AH2 based on sequence and residue number
                    # AH2 is hydrophobic
                    residue_number = np.int64(ndx / self.beads_per_residue)
                    if self.sequence[residue_number] in self.hydrophobic_dict:
                        snap.particles.typeid[idx] = self.getTypeByName['AH2']
                    else:
                        snap.particles.typeid[idx] = self.getTypeByName['AH1']

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
        print(f"AH information (helix copolymer model)")
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
        print(f"  B self hydrophobic                = {self.Bself_hydrophobic}")
        print(f"  B self other                      = {self.Bself_other}")
        print(f"AH R (sigma)                        = {self.r0}")
        print(f"AH RC (sigma)                       = {self.rc}")
        print(f"AH r (sigma)                        = {self.rbond}")
        print(f"AH kbond (kBT/sigma^2)              = {self.kbond}")
        print(f"AH kbend (kBT/rad^2)                = {self.kbend}")
        print(f"AH thetabend (rad)                  = {self.thetabend}")
        print(f"--Helix parameters--")
        print(f"Initial arrangement                 = {self.initial_arrangement}")
        print(f"Sequence                            = {self.sequence}")
        print(f"Effective number of amino acids     = {self.naa}")
        print(f"Number of beads per amino acid      = {self.beads_per_residue}")
        print(f"Helix radius(sigma)                 = {self.radius}")
        print(f"Helix pitch (sigma)                 = {self.pitch}")
        print(f"Helix number of beads               = {self.nbeads}")
        print(f"Helix total length (sigma)          = {self.total_length}")
