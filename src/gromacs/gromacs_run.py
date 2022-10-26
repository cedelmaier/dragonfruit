#!/usr/bin/env python3

# XXX: Put a license here

"""Class for a gromacs run"""

import gc
import logging
import pickle
import os
import sys
import time
import yaml

# Import MDAnalysis as mda
import MDAnalysis as mda

# Numpy and matplotlib stuff
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src', 'lib'))
from stylelib.common_styles import *
from gromacs_seed import GromacsSeed
from gromacs_sim import GromacsSim
from run_base import RunBase

from common import create_datadir
from seed_graph_funcs import *
from sim_graph_funcs import *

class GromacsRun(RunBase):
    def __init__(self, opts, simdir_path = "simulations", simType = GromacsSim, seedType = GromacsSeed):
        if opts.verbose: print(f"GromacsRun::__init__")
        RunBase.__init__(self, opts, simdir_path, simType = GromacsSim, seedType = GromacsSeed)

        # Try to set up the logging facilities correctly
        mda_logger = logging.getLogger('MDAnalysis')
        if mda_logger.hasHandlers():
            mda_logger.handlers = []
        mda.start_logging()
        self.logger = logging.getLogger('MDAnalysis')

        self.verbose = opts.verbose

        # Here is some of the uniqueness for the Gromacs version of this run, this is hardcoded!
        print(f"WARNING: Hardcoded section of GromacsRun looking for simulations/seeds")
        print(f"WARNING: Hardcoded section of GromacsRun looking for simulations/seeds")
        print(f"WARNING: Hardcoded section of GromacsRun looking for simulations/seeds")

        # Code in where to find everything
        self.main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data/aggregated_symlinks/')

        # There are 4 simulations under different conditions, grab them and note what they are done at
        self.sims = {}
        self.sims['polarmonomer_aglipid_rotxN_50mMKCl'] = ('C-cap', 50)
        self.sims['polarmonomer_aglipid_rotxN_150mMKCl'] = ('C-cap', 150)
        self.sims['neutralmonomer_aglipid_rotxN_50mMKCl'] = ('N-cap', 50)
        self.sims['neutralmonomer_aglipid_rotxN_150mMKCl'] = ('N-cap', 150)

        self.named_graphs = {}
        self.named_graphs['zpos']          = r"Z ($\AA$)"
        self.named_graphs['helix']         = r"Helicity (AU)"
        self.named_graphs['tilt']          = r"$\Theta_{u}$ (deg)"
        self.named_graphs['pdip']          = r"$\Theta_{p}$ (deg)"
        self.named_graphs['pmom']          = r"Helix electric dipole magnitude"
        self.named_graphs['hp']            = r"$\Phi_{up}$ (deg)"
        self.named_graphs['zforce']        = r"Force (kJ mol$^{-1}$ nm${^-1}$)"
        self.named_graphs['perptorque']    = r"Helix perpendicular torque (kJ mol$^{-1}$)"

        # Load up the ylow ad yhi dicts
        self.ylow_dict = {"zpos": -30.0,
                          "helix": 0.0,
                          "tilt": 0.0,
                          "pdip": 0.0,
                          "pmom": 0.0,
                          "hp": 0.0,
                          "zforce": -8000.0,
                          "perptorque": -800000.0}
        self.yhi_dict = {"zpos": 30.0,
                          "helix": 1.05,
                          "tilt": 180.0,
                          "pdip": 180.0,
                          "pmom": 40.0,
                          "hp": 180.0,
                          "zforce": 8000.0,
                          "perptorque": 800000.0}

        # Going to load the dataframes with the data too
        self.dfs = {}
        for key,_ in self.sims.items():
            print(f"Loading final data: {key}")
            file_path_current = os.path.join(self.main_path, key, "data", "{}.h5".format(key))
            self.dfs[key] = pd.read_hdf(file_path_current)
            print(self.dfs[key])

        if opts.verbose: print(f"GromacsRun::__init__ return")

    def GraphRun(self):
        r""" Graph this run
        """
        if self.verbose: print(f"GromacsRun::GraphRun")

        self.GraphDistributions()
        self.GraphFinalScatter()

        if self.verbose: print(f"GromacsRun::GraphRun return")

    def GraphDistributions(self):
        r""" Aggregate similar condition simulations
        """
        if self.verbose: print(f"GromacsRun::GraphDistributions")

        colors = mpl.cm.rainbow(np.linspace(0,1,len(self.sims)))
        self.datadir = create_datadir(self.main_path, "data")

        # Loop over graphs
        for key,val in self.named_graphs.items():
            print(f"Distribution: {key}")
           
            fig_dist, ax_dist   = plt.subplots(1, 1, figsize = (16, 9))
            plt.figure(fig_dist)
            # Loop over simulations
            fig_range = (self.ylow_dict[key], self.yhi_dict[key])
            legend_names = []
            xdatas = []
            for sim, simname in self.sims.items():
                xdata = self.dfs[sim][[key]].to_numpy()
                xdatas.append(xdata.flatten())

                legend_names.append("{} {} mM".format(simname[0], simname[1]))

            # Do this as a single histogram for the pretty bar plots
            xdatas = np.array(xdatas)
            xdatas = xdatas.T
            his_ret = ax_dist.hist(xdatas, histtype='bar', color = colors, range = fig_range)
            print(his_ret)
            ax_dist.legend(legend_names)
            ax_dist.set_xlabel(self.named_graphs[key])
            ax_dist.set_ylabel("Count")
            fig_dist.tight_layout()
            plt.savefig("{}/distribution_{}.pdf".format(self.datadir, key), dpi = fig_dist.dpi)

            plt.close()
            gc.collect()

        if self.verbose: print(f"GromacsRun::GraphDistributions return")

    def GraphFinalScatter(self):
        r""" Graph scatter plots of every final variable versus every other final variable
        """
        if self.verbose: print(f"GromacsRun::GraphFinalScatter")

        colors = mpl.cm.rainbow(np.linspace(0,1,len(self.sims)))
        self.combodir = create_datadir(self.datadir, "combinations")

        # Generate all unique pairs of keys
        combinations = list(itertools.combinations(self.named_graphs.keys(), 2))

        # Now loop over them
        for combo in combinations:
            print(f"Combination: {combo}")
            fig_scatter, ax_scatter = plt.subplots(1, 1, figsize = (16, 9))

            key1 = combo[0]
            key2 = combo[1]

            legend_names = []
            color_count = 0
            for sim, simname in self.sims.items():
                xdata = self.dfs[sim][[key1]].to_numpy().flatten()
                ydata = self.dfs[sim][[key2]].to_numpy().flatten()

                ax_scatter.scatter(xdata, ydata, color = colors[color_count])
                legend_names.append("{} {} mM".format(simname[0], simname[1]))

                color_count += 1

            ax_scatter.set_xlim([self.ylow_dict[key1], self.yhi_dict[key1]])
            ax_scatter.set_ylim([self.ylow_dict[key2], self.yhi_dict[key2]])
            ax_scatter.legend(legend_names)
            ax_scatter.set_xlabel(self.named_graphs[key1])
            ax_scatter.set_ylabel(self.named_graphs[key2])
            fig_scatter.tight_layout()
            plt.savefig("{}/combination_{}_{}.pdf".format(self.combodir, key1, key2), dpi = fig_scatter.dpi)

            # clean up
            plt.close(fig_scatter)
            gc.collect()

        if self.verbose: print(f"GromacsRun::GraphFinalScatter return")