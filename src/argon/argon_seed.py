# XXX: Put a license here

""" Class for a single argon seed """

import os
import random
import sys
import yaml

import gsd.hoomd
import numpy as np
import pandas as pd

# Analysis imports
import freud
from freud import box

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from seed_base import SeedBase

class ArgonSeed(SeedBase):
    def __init__(self, path, opts):
        print(f"ArgonSeed::__init__")

        # Run the base class initialization
        SeedBase.__init__(self, path, opts)

        # Now we do any other kind of initialization we require
        self.ReadData()

        # Set the RNG state (VERY IMPORTANT)
        random.seed(self.nseed)

        print(f"ArgonSeed::__init__ return")

    def ReadData(self):
        r""" Read the data from a YAML file for a single seed
        """
        print(f"ArgonSeed::ReadData")
        # XXX: Figure out what all needs to be read in here!
        self.kT             = np.float64(self.default_yaml['simulation']['kT'])
        self.deltatau       = np.float64(self.default_yaml['simulation']['deltatau'])
        self.nseed          = np.int32(np.float64(self.default_yaml['simulation']['seed']))
        self.compute_mode   = self.default_yaml['simulation']['mode']

        print(f"ArgonSeed::ReadData return")

    def PrintInformation(self, snap):
        r""" Print parameter information for ArgonSeed
        """
        # XXX: Figure out how to print out the informatio here (see septin_seed.py)

    def Configure(self, snap):
        r""" Configure the simulation snapshot for HOOMD
        """
        print(f"ArgonSeed::Configure")

        # XXX: Figure out how to configure the system, see septin_seed.py

        print(f"ArgonSeed::Configure return")
