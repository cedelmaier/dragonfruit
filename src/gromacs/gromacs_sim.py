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

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from stylelib.common_styles import septin_runs_stl
from gromacs_seed import GromacsSeed
from sim_base import SimulationBase

from seed_graph_funcs import *

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

        # Define the functions we want to graph on a per-seed basis
        #self.graph_perseed_functions = [ graph_seed_zpos_wheads,
        #                                 graph_seed_helicity,
        #                                 graph_seed_globaltilt,
        #                                 graph_seed_pdipoletilt,
        #                                 graph_seed_helixpdipoleangle ]

        self.graph_perseed_trajectory   = [ graph_seed_zpos_wheads,
                                            graph_seed_helicity ]
        self.graph_perseed_tilts        = [ graph_seed_globaltilt,
                                            graph_seed_pdipoletilt,
                                            graph_seed_helixpdipoleangle ]
        self.graph_perseed_forces       = [ graph_seed_zforce,
                                            graph_seed_perptorque ]

        # Combine the above into a single thing we can graph and have text on and what have you
        self.graph_groups = {}
        self.graph_groups["trajectory"] = self.graph_perseed_trajectory
        self.graph_groups["angles"] = self.graph_perseed_tilts
        self.graph_groups["forces"] = self.graph_perseed_forces


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
        if self.verbose: print(f"GromacsSim::GraphSimulation")
        # A couple of quick checks to make sure data directories are setup correctly
        if not self.opts.datadir: self.opts.datadir = create_datadir(self.opts.workdir)
        if self.opts.workdir != self.sim_path:
            self.sim_datadir = create_datadir(self.opts.datadir, datadir_name = "{}_data".format(self.name))
        else:
            self.sim_datadir = self.opts.datadir

        #self.GraphDynamicData()
        self.GraphGroups()

        if self.verbose: print(f"GromacsSim::GraphSimulation return")

    def GraphGroups(self):
        r""" Graph groups of data for simulation
        """
        for key,val in self.graph_groups.items():
            n_graphs = len(val)
            fig, axarr = plt.subplots(n_graphs, 2, figsize=(25,16))
            colors = mpl.cm.rainbow(np.linspace(0,1,len(self.seeds)))

            # Get the timing information
            timepoints_traj = self.seeds[-1].master_time_df.index/1000.0
            timepoints_forces = self.seeds[-1].master_forces_df.index/1000.0

            # Zip togther the list like we had before to use the colors
            for graph, axr in zip(val, axarr):
                # Get the average data
                yarr_avg = None
                num_seeds = 0 # Used later if successful
                # Create a graph for each seed color
                for sd, col in zip(self.seeds, colors):
                    [ylow, yhi, yarr] = graph(sd, axr[0], color=col, xlabel = False)
                    if yarr_avg is None:
                        yarr_avg = np.zeros_like(yarr)
                    yarr_avg = np.add(yarr_avg, yarr)
                    num_seeds += 1

                # Actually do the average
                yarr_avg /= np.float32(num_seeds)

                # Check for different behavior based on the number of returns
                if key == "trajectory" or key == "angles":
                    timepoints = timepoints_traj
                else:
                    timepoints = timepoints_forces

                if yarr_avg.ndim == 2:
                    axr[1].plot(timepoints, yarr_avg[0][:], color = "slategrey")
                    axr[1].plot(timepoints, yarr_avg[1][:], color = "slategrey")
                    axr[1].plot(timepoints, yarr_avg[2][:])
                else:
                    axr[1].plot(timepoints, yarr_avg)
                axr[1].set_ylim([ylow, yhi])

            # Set axis labels
            axarr[-1,1].set_xlabel(r'Time (ns)')
            axarr[-1,0].set_xlabel(r'Time (ns)')

            fig.tight_layout()
            plt.savefig("{}_{}.pdf".format(os.path.join(self.sim_datadir, key), self.name), dpi=fig.dpi)

            # Clean up
            fig.clf()
            plt.close()
            gc.collect()
            #mpl.rcdefaults()


    def GraphDynamicData(self):
        r""" Graph dynamic (time) data for simulation
        """
        # For now, setup the 3 different graphs
        fig, axarr = plt.subplots(2, 2, figsize=(25,16))
        colors = mpl.cm.rainbow(np.linspace(0,1,len(self.seeds)))

        # Get the timing information (in ns)
        timepoints = self.seeds[-1].master_time_df.index/1000.0

        # Zip together the list like we had before to be clever
        # The first loop zips together the functions we want to call with the axis array we generated
        for graph, axr in zip(self.graph_perseed_trajectory, axarr):
            # Get the average array all setup correctly
            yarr_avg = None
            num_seeds = 0 # Use later for successful determination
            # The second loop creates different colors for each seed
            for sd, col in zip(self.seeds, colors):
                [ylow, yhi, yarr] = graph(sd, axr[0], color=col, xlabel = False) # Actal magic of the grahping call
                if yarr_avg is None:
                    yarr_avg = np.zeros_like(yarr)
                yarr_avg = np.add(yarr_avg, yarr)
                num_seeds += 1

            yarr_avg /= np.float32(num_seeds)
            # Check for different behavior if we have multiple arrays
            if yarr_avg.ndim == 2:
                axr[1].plot(timepoints, yarr_avg[0][:], color = "slategrey")
                axr[1].plot(timepoints, yarr_avg[1][:], color = "slategrey")
                axr[1].plot(timepoints, yarr_avg[2][:])
            else:
                axr[1].plot(timepoints, yarr_avg)
            axr[1].set_ylim([ylow, yhi])

        # XXX: Legend isn't plotting properly figure this out some other time
        #axarr[1,0].legend(loc='center right', bbox_to_anchor=(1.0, 0.0))

        axarr[-1,1].set_xlabel(r'Time (ns)')
        axarr[-1,0].set_xlabel(r'Time (ns)')

        fig.tight_layout()
        plt.savefig("{}_{}.pdf".format(os.path.join(self.sim_datadir, 'dynamicdata'), self.name), dpi=fig.dpi)

        # Clean up
        fig.clf()
        plt.close()
        gc.collect()
        mpl.rcdefaults()

