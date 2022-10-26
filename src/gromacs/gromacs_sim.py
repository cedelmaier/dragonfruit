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

from common import create_datadir
from seed_graph_funcs import *
from sim_graph_funcs import *

# Move thit into a subclass at some point
def shaded_error(ax, x, y, error, alpha, color, linestyle = 'solid', label = None):
    r""" Plot a shaded error bar with colors
    """
    ax.plot(x, y, color = color, linestyle = linestyle, label = label)
    ax.fill_between(x, y-error, y+error,
                    alpha = alpha, edgecolor = color, facecolor = color, linestyle = linestyle)

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
        # Per-seed graphs that we want to make (single graphs for quality)
        self.graph_perseed = {}
        self.graph_perseed['zpos']          = graph_seed_zpos_wheads
        self.graph_perseed['helix']         = graph_seed_helicity
        self.graph_perseed['tilt']          = graph_seed_globaltilt
        self.graph_perseed['pdip']          = graph_seed_pdipoletilt
        self.graph_perseed['pmom']          = graph_seed_pdipolemoment
        self.graph_perseed['hp']            = graph_seed_helixpdipoleangle
        self.graph_perseed['zforce']        = graph_seed_zforce
        self.graph_perseed['perptorque']    = graph_seed_perptorque

        self.graph_yaxis_title = {}
        self.graph_yaxis_title['zpos']          = r"Z ($\AA$)"
        self.graph_yaxis_title['helix']         = r"Helicity (AU)"
        #self.graph_yaxis_title['tilt']          = r"Helix tilt (deg)"
        self.graph_yaxis_title['tilt']          = r"$\Theta_{u}$ (deg)"
        #self.graph_yaxis_title['pdip']          = r"Helix electric dipole tilt (deg)"
        self.graph_yaxis_title['pdip']          = r"$\Theta_{p}$ (deg)"
        self.graph_yaxis_title['pmom']          = r"Helix electric dipole magnitude"
        #self.graph_yaxis_title['hp']            = r"Helix-Helix electric dipole angle (deg)"
        self.graph_yaxis_title['hp']            = r"$\Phi_{up}$ (deg)"
        self.graph_yaxis_title['zforce']        = r"Force (kJ mol$^{-1}$ nm${^-1})"
        self.graph_yaxis_title['perptorque']    = r"Helix perpendicular torque (kJ mol$^{-1}$)"

        # Set up the final data information
        self.final_data = {}
        self.ylow_dict = {}
        self.yhi_dict = {}
        for key,_ in self.graph_yaxis_title.items():
            self.final_data[key] = []
        self.dfs = []

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

        # Categorize the final state of the seed by tilt
        self.CategorizeSeeds()

        print(f"---- Simulation {self.name} analyzed ----")

        if self.verbose: print(f"GromacsSim::Analyze return")

    def CategorizeSeeds(self):
        r""" Categorize the seeds by tilt
        """
        if self.verbose: print(f"GromacsSim::CategorizeSeeds")

        print(f"WARNING: Categorize Seeds was done in a hurray, check later!")
        print(f"WARNING: Categorize Seeds was done in a hurray, check later!")
        print(f"WARNING: Categorize Seeds was done in a hurray, check later!")

        tilt_criteria = 50.0
        end_times = 175.0

        for sd in self.seeds:
            print(f"Categorizing seed: {sd.label}")
            global_tilt = sd.master_time_df['global_tilt'].to_numpy().flatten()
            times = sd.master_time_df.index
            times = times/1000.0
            flat_times = times.to_numpy().flatten()

            final_indices = np.where(flat_times > end_times)

            final_tilt_seed = np.mean(global_tilt[final_indices])
            is_tilted = (final_tilt_seed < tilt_criteria) or (final_tilt_seed > 180.0 - tilt_criteria)
            print(f"  Criteria: {is_tilted}, Final tilt: {final_tilt_seed}")

        if self.verbose: print(f"GromacsSim::CategorizeSeeds return")

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

        self.GraphPerseed()
        self.GraphGroups()
        self.GraphDistributions()
        self.GraphFinalScatter()

        # XXX Move this into a separate write data section for the future
        self.final_df = pd.concat(self.dfs, axis=1)
        print(self.final_df)

        if self.verbose: print(f"GromacsSim::GraphSimulation return")

    def GraphPerseed(self):
        r""" Graph perseed information in a single high quality plot
        """
        if self.verbose: print(f"GromacsSim::GraphPerseed")

        individual_dir = create_datadir(self.sim_datadir, "individual")
        self.individual_dir = individual_dir

        # Loop over the seeds and grab the information
        nseeds = len(self.seeds)
        print(f"Graphing individual quantities ({individual_dir})")
        for key,val in self.graph_perseed.items():
            print(f"  Graph: {key}")
            fig, axarr                              = plt.subplots(1, 2, figsize = (16,9))
            fig_seeds, axarr_seeds                  = plt.subplots(1, 1, figsize = (16, 9))
            fig_mean, axarr_mean                    = plt.subplots(1, 1, figsize = (16, 9))
            graph = val # Alias graph
            colors = mpl.cm.rainbow(np.linspace(0,1,nseeds))

            # Get the timing information
            timepoints_traj = self.seeds[-1].master_time_df.index/1000.0
            timepoints_forces = self.seeds[-1].master_forces_df.index/1000.0

            # different stuff has different plotting requirements
            if key == "zforce" or key == "perptorque":
                timepoints = timepoints_forces
                dolegend = False
            else:
                timepoints = timepoints_traj
                dolegend = True

            # Get the average data set up
            yarr_list = []
            ylow_list = []
            yhi_list = []
            for sd, col in zip(self.seeds, colors):
                [ylow, yhi, yarr] = graph(sd, axarr[0], color=col, xlabel = False, dolegend = dolegend)
                # Repeat the graph on the seeds only version
                graph(sd, axarr_seeds, color = col, xlabel = False, dolegend = dolegend)
                yarr_list.append(yarr)
                ylow_list.append(ylow)
                yhi_list.append(yhi)

                self.ylow_dict[key] = ylow
                self.yhi_dict[key] = yhi

                # Grab the information on the single seed
                # XXX: Figure out how to make this easier, or move the data harvesting to an analysis location
                yarr_nparray = np.array(yarr)
                if yarr_nparray.ndim == 2:
                    if key == "zforce" or key == "perptorque":
                        self.final_data[key].append(yarr_nparray[:,-1])
                    else:
                        self.final_data[key].append(yarr_nparray[:,-1])
                else:
                    self.final_data[key].append(yarr_nparray[-1])

            # Get the mean and STD
            y_mean = np.mean(np.array(yarr_list), axis=0)
            y_std = np.std(np.array(yarr_list), axis=0, ddof=1)

            # Check the plot for things like upper/lower leaflet stuff
            if y_mean.ndim == 2:
                # Check if we are doing forces/torques, or just the membrane
                if key == "zforce" or key == "perptorque":
                    shaded_error(axarr[1], timepoints, y_mean[0][:], y_std[0][:], alpha = 0.5, color = 'b', linestyle = 'solid', label = 'All')
                    shaded_error(axarr[1], timepoints, y_mean[1][:], y_std[1][:], alpha = 0.5, color = 'c', linestyle = 'dotted', label = 'Electrostatic')
                    shaded_error(axarr[1], timepoints, y_mean[2][:], y_std[2][:], alpha = 0.5, color = 'm', linestyle = 'dashed', label = 'Non-electrostatic')
                    axarr[1].legend()

                    # Do the average only plot too
                    shaded_error(axarr_mean, timepoints, y_mean[0][:], y_std[0][:], alpha = 0.5, color = 'b', linestyle = 'solid', label = 'All')
                    shaded_error(axarr_mean, timepoints, y_mean[1][:], y_std[1][:], alpha = 0.5, color = 'c', linestyle = 'dotted', label = 'Electrostatic')
                    shaded_error(axarr_mean, timepoints, y_mean[2][:], y_std[2][:], alpha = 0.5, color = 'm', linestyle = 'dashed', label = 'Non-electrostatic')
                    axarr_seeds.legend()
                else:
                    shaded_error(axarr[1], timepoints, y_mean[0][:], y_std[0][:], alpha = 0.5, color = 'slategrey')
                    shaded_error(axarr[1], timepoints, y_mean[1][:], y_std[1][:], alpha = 0.5, color = 'slategrey')
                    shaded_error(axarr[1], timepoints, y_mean[2][:], y_std[2][:], alpha = 0.5, color = 'b')

                    shaded_error(axarr_mean, timepoints, y_mean[0][:], y_std[0][:], alpha = 0.5, color = 'slategrey')
                    shaded_error(axarr_mean, timepoints, y_mean[1][:], y_std[1][:], alpha = 0.5, color = 'slategrey')
                    shaded_error(axarr_mean, timepoints, y_mean[2][:], y_std[2][:], alpha = 0.5, color = 'b')
            else:
                shaded_error(axarr[1], timepoints, y_mean, y_std, alpha = 0.5, color = 'b')
                shaded_error(axarr_mean, timepoints, y_mean, y_std, alpha = 0.5, color = 'b')

            axarr[1].set_ylim([ylow, yhi])
            axarr_mean.set_ylim([ylow, yhi])

            # Set axis labels
            axarr[0].set_xlabel(r'Time (ns)')
            axarr[1].set_xlabel(r'Time (ns)')

            axarr_seeds.set_xlabel(r'Time(ns)')
            axarr_mean.set_xlabel(r'Time(ns)')

            axarr_mean.set_ylabel(self.graph_yaxis_title[key])

            # Save the figure
            plt.figure(fig)
            fig.tight_layout()
            plt.savefig("{}/{}_{}.pdf".format(individual_dir, key, self.name), dpi = fig.dpi)
            fig.clf()

            plt.figure(fig_seeds)
            fig_seeds.tight_layout()
            plt.savefig("{}/{}_{}_seeds.pdf".format(individual_dir, key, self.name), dpi = fig_seeds.dpi)
            fig_seeds.clf()

            plt.figure(fig_mean)
            fig_mean.tight_layout()
            plt.savefig("{}/{}_{}_mean.pdf".format(individual_dir, key, self.name), dpi = fig_mean.dpi)
            fig_mean.clf()

            plt.close(fig)
            plt.close(fig_seeds)
            plt.close(fig_mean)
            gc.collect()

        if self.verbose: print(f"GromacsSim::GraphPerseed return")

    def GraphGroups(self):
        r""" Graph groups of data for simulation
        """
        if self.verbose: print(f"GromacsSim::GraphGroups")
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
        if self.verbose: print(f"GromacsSim::GraphGroups return")

    def GraphDistributions(self):
        r""" Graph distributions of variables
        """
        if self.verbose: print(f"GromacsSim::GraphDistributions")

        for key,_ in self.graph_perseed.items():
            print(f"  Distribution: {key}")
            fig_dist, ax_dist   = plt.subplots(1, 1, figsize = (16, 9))

            if key == "zforce" or key == "perptorque":
                xdata = np.array(self.final_data[key])[:,0]
            elif key == "zpos":
                xdata = np.array(self.final_data[key])[:,2]
            else:
                xdata = np.array(self.final_data[key])

            plt.figure(fig_dist)
            fig_range = (self.ylow_dict[key], self.yhi_dict[key])
            (hist, bin_mids, npoints) = graph_sim_histogram(xdata, ax_dist, mtitle = key, ytitle = "Count", xtitle = self.graph_yaxis_title[key], xrange = fig_range)
            fig_dist.tight_layout()
            plt.savefig("{}/{}_distribution_{}.pdf".format(self.individual_dir, key, self.name), dpi = fig_dist.dpi)

            # Save off the data!
            df = pd.DataFrame({key: xdata})
            self.dfs.append(df)

            # Clean up
            plt.close(fig_dist)
            gc.collect()

        if self.verbose: print(f"GromacsSim::GraphDistributions return")

    def GraphFinalScatter(self):
        r""" Graph scatter plots of every final variable versus every other final variable
        """
        if self.verbose: print(f"GromacsSim::GraphFinalScatter")

        # Create a new directory for the results
        combodir = create_datadir(self.sim_datadir, "combinations")

        # Generate all unique pairs of keys
        combinations = list(itertools.combinations(self.graph_perseed.keys(), 2))

        # Now loop over them
        for combo in combinations:
            print(f"  Combination: {combo[0]}, {combo[1]}")
            fig_scatter, ax_scatter = plt.subplots(1, 1, figsize=(16,9))

            key1 = combo[0]
            key2 = combo[1]

            if key1 == "zforce" or key1 == "perptorque":
                xdata = np.array(self.final_data[key1])[:,0]
            elif key1 == "zpos":
                xdata = np.array(self.final_data[key1])[:,2]
            else:
                xdata = np.array(self.final_data[key1])

            if key2 == "zforce" or key2 == "perptorque":
                ydata = np.array(self.final_data[key2])[:,0]
            elif key2 == "zpos":
                ydata = np.array(self.final_data[key2])[:,2]
            else:
                ydata = np.array(self.final_data[key2])

            print(xdata)
            print(ydata)

            plt.figure(fig_scatter)
            graph_sim_scatter(xdata, ydata, ax_scatter, mtitle = "{}_{}".format(key1, key2), xtitle = self.graph_yaxis_title[key1], ytitle = self.graph_yaxis_title[key2])
            ax_scatter.set_xlim([self.ylow_dict[key1], self.yhi_dict[key1]])
            ax_scatter.set_ylim([self.ylow_dict[key2], self.yhi_dict[key2]])
            fig_scatter.tight_layout()
            plt.savefig("{}/{}_{}_{}.pdf".format(combodir, key1, key2, self.name), dpi = fig_scatter.dpi)

            # clean up
            plt.close(fig_scatter)
            gc.collect()

        if self.verbose: print(f"GromacsSim::GraphFinalScatter return")
