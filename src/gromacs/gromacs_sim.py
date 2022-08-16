#!/usr/bin/env python3

# XXX: Put a license here

"""Class for a single gromacs simulation"""

import gc
import logging
import pickle
import os
import sys
import time
import yaml

# Import MDAnalysis as mda
import MDAnalysis as mda

# MDTraj has a native DSSP module
import mdtraj as md

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from stylelib.common_styles import septin_runs_stl
from gromacs_seed import GromacsSeed
from sim_base import SimulationBase

class GromacsSim(SimulationBase):
    def __init__(self, path, opts, seedType=GromacsSeed):
        if opts.verbose: print(f"GromacsSim::__init__")
        SimulationBase.__init__(self, path, opts, seedType=seedType)
        # Try to set up the logging facilities correctly
        mda_logger = logging.getLogger('MDAnalysis')
        if mda_logger.hasHandlers():
            mda_logger.handlers = []
        mda.start_logging()
        self.logger = logging.getLogger('MDAnalysis')

        self.verbose = opts.verbose
        if self.verbose: print("GromacsSim::__init__ return")

    def Analyze(self):
        r"""Analyze all of the underlying seeds
        """
        if self.verbose: print(f"GromacsSim::Analyze")
        for sd in self.seeds: sd.Analyze()
        print(f"---- Simulation {self.name} analyzed ----")
        if self.verbose: print(f"GromacsSim::Analyze return")

    def GraphSimulation(self):
        r""" Graph simulation parameters and variables
        """
        # A couple of quick checks to make sure data directories are setup correctly
        if not self.opts.datadir: self.opts.datadir = create_datadir(self.opts.workdir)
        if self.opts.workdir != self.sim_path:
            self.sim_datadir = create_datadir(self.opts.datadir, datadir_name = "{}_data".format(self.name))
        else:
            self.sim_datadir = self.opts.datadir

        self.GraphDynamicData()

    def GraphDynamicData(self):
        r""" Graph dynamic (time) data for simulation
        """
        # For now, setup the 3 different graphs

