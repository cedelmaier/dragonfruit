# XXX: Copyright info here

"""Class for block-style AH domains"""

import hoomd
import hoomd.md as md
import gsd.hoomd

import enum
import itertools
import os
import random
import sys
import yaml

import numpy as np
from ahdomain import ahdomain

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from common import *

# Binding state enum
class BindingState(enum.Enum):
    init = -1
    free = 0
    near = 1
    surface = 2
    intermediate = 3
    deep = 4

class SingleAHDomain(object):
    r""" Single Block AH Domain class for analysis
    """

    def __init__(self, idx, nbeads):
        self.idx = idx                  # Which AH domain we are
        self.idx_start = -1             # Starting bead for this AH
        self.nbeads = nbeads            # Number of beads in the chain

        # State information
        self.current_state = BindingState.free  # Start in the free binding state
                                        # 1 for near surface, 2 for surface interact,
                                        # 3 for intermediate, 4 for deep
        self.current_state_counter = 1  # Counter to keep track of how many frames we've been in the current state
        self.lifetime_dict = {}         # Dictionary to keep track of global lifetimes
        self.lifetime_dict[BindingState.free] = []
        self.lifetime_dict[BindingState.near] = []
        self.lifetime_dict[BindingState.surface] = []
        self.lifetime_dict[BindingState.intermediate] = []
        self.lifetime_dict[BindingState.deep] = []
        self.lifetime_subdivision_dict = {}

        # Indexes for various things
        self.index_all = None
        self.index_ah1 = None
        self.index_ah2 = None

    def UpdateBindingState(self, current_binding_state, current_time_division):
        # Unfortunately have to do some fancy work for the subdivision information
        if current_time_division not in self.lifetime_subdivision_dict:
            self.lifetime_subdivision_dict[current_time_division] = {}
            self.lifetime_subdivision_dict[current_time_division][BindingState.free] = []
            self.lifetime_subdivision_dict[current_time_division][BindingState.near] = []
            self.lifetime_subdivision_dict[current_time_division][BindingState.surface] = []
            self.lifetime_subdivision_dict[current_time_division][BindingState.intermediate] = []
            self.lifetime_subdivision_dict[current_time_division][BindingState.deep] = []

        # Check to see if there is a state change
        if self.current_state != current_binding_state:
            # Use reflection in the enum to figure out how to do this
            self.lifetime_dict[self.current_state].append(self.current_state_counter)
            self.lifetime_subdivision_dict[current_time_division][self.current_state].append(self.current_state_counter)

            # Update the state and reset the counter
            self.current_state = current_binding_state
            self.current_state_counter = 1
        else:
            # Increment the frame counter
            self.current_state_counter += 1


    def PrintInformation(self):
        print(f"Idx: {self.idx}")
        print(f"  idx_start: {self.idx_start}")
        print(f"  nbeads: {self.nbeads}")
        print(f"  current_state: {self.current_state}")
        print(f"  index_all: {self.index_all}")
        print(f"  index_ah1: {self.index_ah1}")
        print(f"  index_ah2: {self.index_ah2}")
        print(f"    free_lifetime:          {self.lifetime_dict[BindingState.free]}")
        print(f"    near_lifetime:          {self.lifetime_dict[BindingState.near]}")
        print(f"    surface_lifetime:       {self.lifetime_dict[BindingState.surface]}")
        print(f"    intermediate_lifetime:  {self.lifetime_dict[BindingState.intermediate]}")
        print(f"    deep_lifetime:          {self.lifetime_dict[BindingState.deep]}")

class BlockAHDomain(ahdomain):
    def __init__(self, bead_size, yaml_file):
        super().__init__(bead_size, yaml_file)

        # Set constants we will need for this simulation
        self.nbeads           = np.int32(np.float64(yaml_file['ah_domain']['nbeads']))
        self.nrepeat          = np.int32(np.float64(yaml_file['ah_domain']['nrepeat']))
        self.thetabend        = np.float64(yaml_file['ah_domain']['thetabend'])

        # Set constants we will need
        self.bead_size = 0.75*bead_size
        self.mass_per_bead = self.mass / self.nbeads
        self.r0 = self.bead_size
        self.rc = 2.0*self.r0
        self.rbond = self.r0

        # Initial configuration information for analysis
        self.analysis_init = False
        self.getTypeIDbyName = {}
        self.getNamebyTypeID = {}

    # Create some number of block copolyerms in the simulation
    def CreateAH(self, snap):
        # Early bail if no AH domains
        if self.nah == 0:
            return

        # We have additional types that we must setup before creating particles
        p_ntypes = len(snap.particles.types)
        snap.particles.types = snap.particles.types + ['AH1', 'AH2']
        self.getTypeIDbyName = {'AH1': p_ntypes+0, 'AH2': p_ntypes+1}
        b_ntypes = len(snap.bonds.types)
        snap.bonds.types = snap.bonds.types + ['ahbond']
        self.getTypeIDbyName['ahbond'] = b_ntypes
        a_ntypes = len(snap.angles.types)
        snap.angles.types = snap.angles.types + ['ahbend']
        self.getTypeIDbyName['ahbend'] = a_ntypes

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
                        snap.particles.typeid[idx] = self.getTypeIDbyName['AH1']
                    else:
                        snap.particles.typeid[idx] = self.getTypeIDbyName['AH2']
                    snap.particles.mass[idx] = self.mass_per_bead

                    ndx += 1

                # Now set up the bond information
                n_newbonds = self.nbeads - 1
                ah_start_bdx = snap.bonds.N
                bdx = ah_start_bdx
                snap.bonds.N = snap.bonds.N + n_newbonds
                for idx in range(ah_start_nx + self.nbeads*ahdx, ah_start_nx + self.nbeads*ahdx + self.nbeads-1):
                    snap.bonds.typeid[bdx] = self.getTypeIDbyName['ahbond']
                    snap.bonds.group[bdx] = [idx, idx+1]
                    bdx += 1

                # Now set up the angle information
                n_newangles = self.nbeads - 2
                ah_start_adx = snap.angles.N
                adx = ah_start_adx
                snap.angles.N = snap.angles.N + n_newangles
                for idx in range(ah_start_nx + self.nbeads*ahdx, ah_start_nx + self.nbeads*ahdx + self.nbeads - 2):
                    snap.angles.typeid[adx] = self.getTypeIDbyName['ahbend']
                    snap.angles.group[adx] = [idx, idx+1, idx+2]
                    adx += 1

                # Bail for the while loop
                if overlap == 0:
                    break
       
    def PrintInformation(self, snap):
        # Print out AH information
        print(f"--------")
        print(f"AH information (block copolymer model)")
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
        print(f"  B self (hydrophobic)              = {self.Bself_hydrophobic}")
        print(f"  B self (other)                    = {self.Bself_other}")
        print(f"AH R (sigma)                        = {self.r0}")
        print(f"AH RC (sigma)                       = {self.rc}")
        print(f"AH r (sigma)                        = {self.rbond}")
        print(f"AH kbond (kBT/sigma^2)              = {self.kbond}")
        print(f"AH kbend (kBT/rad^2)                = {self.kbend}")
        print(f"AH thetabend (rad)                  = {self.thetabend}")

    def ConfigureAnalysis(self, snap):
        r"""Configure the analysis for the trajectory
        """

        if self.analysis_init:
            print(f"BlockAHDomain: Analysis already set up!")
            return

        # This is for the copolymer model, some linear number of bonds
        self.ah_start_nx = snap.particles.N - (self.nbeads * self.nah)

        # Get the TypeID by name or by type
        for idx,val in enumerate(snap.particles.types):
            self.getTypeIDbyName[val] = idx
            self.getNamebyTypeID[idx] = val

        # Loop over number of AH domains
        self.ahdomains = [SingleAHDomain(ahdx, self.nbeads) for ahdx in range(self.nah)]
        for ahdx in range(self.nah):
            subahdomain = self.ahdomains[ahdx]
            subahdomain.idx_start = self.ah_start_nx + self.nbeads*ahdx

            # Assign particles IDs to the bead list
            # XXX There is probably a better way to do this is the correct boolean indexing and list operations
            subahdomain.index_all = np.arange(self.ah_start_nx + self.nbeads*ahdx, self.ah_start_nx + self.nbeads*(ahdx+1))
            subahdomain.index_ah1 = subahdomain.idx_start + \
                    np.where(snap.particles.typeid[subahdomain.index_all] == self.getTypeIDbyName['AH1'])[0]
            subahdomain.index_ah2 = subahdomain.idx_start + \
                    np.where(snap.particles.typeid[subahdomain.index_all] == self.getTypeIDbyName['AH2'])[0]

        self.analysis_init = True

    def GetSubAHDomain(self, ahdx):
        r""" Get the AH domain (subdomain) for this index
        """
        return self.ahdomains[ahdx]

