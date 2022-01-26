# XXX: Copyright info here

"""Class for Simple Spheres in HOOMD"""

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

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from common import *

class SimpleSpheres(object):
    def __init__(self, bead_size, yaml_file):
        # Simple parameters
        self.nspheres   = np.int32(np.float64(yaml_file['simple_spheres']['nspheres']))
        self.mass       = np.float64(yaml_file['simple_spheres']['mass'])
        self.gamma      = np.float64(yaml_file['simple_spheres']['gamma'])
        self.is_init    = yaml_file['simple_spheres']['is_init']
        self.A          = np.float64(yaml_file['simple_spheres']['A'])
        self.bead_size  = bead_size

        # Make sure we get the sizes right
        self.r0 = self.bead_size
        self.rc = 2.0 * self.r0

        # Initial configuration information for analysis
        self.analysis_init = False
        self.getTypebyName = {}
        self.getNamebyType = {}

    def InitSimpleSpheres(self, snap):
        if self.is_init:
            print(f"Creating simple spheres from file")
            print(f"Not actually sure what to do here, exiting!")
            sys.exit(1)
        else:
            self.CreateSimpleSpheres(snap)
        self.is_init = True

    def CreateSimpleSpheres(self, snap):
        # Early bail if there are no spheres
        if self.nspheres == 0:
            return

        # We have additional types that we should keep track of
        p_ntypes = len(snap.particles.types)
        snap.particles.types = snap.particles.types + ['SS1']
        self.getTypeByName = {'SS1': p_ntypes+0}

        # This is for single spheres, so pretty easy
        spheres_start_nx = snap.particles.N
        snap.particles.N = snap.particles.N + self.nspheres
        rboxxy = snap.configuration.box[0]/2.0
        rboxz = snap.configuration.box[2]/2.0
        # Loop over and insert spheres
        for isphere in range(self.nspheres):
            idx = spheres_start_nx + isphere
            xrand = random.uniform(-1.0*rboxxy, 1.0*rboxxy)
            yrand = random.uniform(-1.0*rboxxy, 1.0*rboxxy)
            zrand = random.uniform(-1.0*rboxz, 1.0*rboxz)

            r1 = np.array([xrand, yrand, zrand])
            snap.particles.position[idx] = r1
            snap.particles.typeid[idx] = self.getTypeByName['SS1']
            snap.particles.mass[idx] = self.mass

    def PrintInformation(self, snap):
        # Print out the simple sphere information
        print(f"--------")
        print(f"Simple Spheres")
        print(f"nspheres                            = {self.nspheres}")
        print(f"mass                                = {self.mass}")
        print(f"gamma                               = {self.gamma}")
