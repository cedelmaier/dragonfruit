# XXX: Copyright info here

"""Base class for run object"""

import ast
import os
import re
import yaml

import numpy as np

from seed_base import SeedBase
from sim_base import SimulationBase

class RunBase(object):
    def __init__(self, opts, simdir_path = "simulations", simType = SimulationBase, seedType = SeedBase):
        self.opts = opts
        if self.opts.verbose: print(f"RunBase::__init__")

        self.run_path = os.path.abspath(opts.workdir)
        self.run_name = self.run_path.split('/')[-1]
        self.simdir_path = os.path.join(self.run_path, simdir_path)

        # Objects to be filled later
        self.sims = []
        self.p_names = []
        self.p_range = []

        self.simType = simType
        self.seedType = seedType

        print(f"WARNING: Run Base is not yet fully implemented, just a hack to get Gromacs to graph correctly")
        print(f"WARNING: Run Base is not yet fully implemented, just a hack to get Gromacs to graph correctly")
        print(f"WARNING: Run Base is not yet fully implemented, just a hack to get Gromacs to graph correctly")

        if self.opts.verbose: print(f"RunBase::__init__ return")
