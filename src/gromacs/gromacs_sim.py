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
        # Add one that is just the alpharmsd
        self.graph_perseed['alpharmsd']     = graph_seed_alpharmsd

        self.named_graphs = {}
        self.named_graphs['zpos']           = r"Z ($\AA$)"
        self.named_graphs['helix']          = r"Helicity (AU)"
        #self.named_graphs['tilt']           = r"$\Theta_{u}$ (deg)"
        self.named_graphs['tilt']           = r"$\Theta_{u}$ (deg)"
        self.named_graphs['pdip']           = r"Tilt (deg)"
        self.named_graphs['pmom']           = r"Helix electric dipole magnitude"
        self.named_graphs['hp']             = r"$\Phi_{up}$ (deg)"
        self.named_graphs['zforce']         = r"Force (kJ mol$^{-1}$ nm${^-1})"
        self.named_graphs['perptorque']     = r"Helix perpendicular torque (kJ mol$^{-1}$)"
        self.named_graphs['alpharmsd']      = r"Helical content (AU)"

        # Set up block averaging information
        self.block_size = 50000.0 # 50 ns blocks
        self.block_average = {}
        self.ylow_dict = {}
        self.yhi_dict = {}

        #self.graph_perseed_trajectory   = [ graph_seed_zpos_wheads,
        #                                    graph_seed_helicity ]
        self.graph_perseed_trajectory   = [ graph_seed_zpos_wheads,
                                            graph_seed_alpharmsd ]
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

        # Look for a stability analysis file for ourselves...
        self.stability_filename = None
        if os.path.isfile(os.path.join(self.opts.workdir,"stability_data.csv")):
            self.stability_filename = os.path.join(self.opts.workdir, "stability_data.csv")

        if self.verbose: print("GromacsSim::__init__ return")

    def Analyze(self):
        r"""Analyze all of the underlying seeds
        """
        if self.verbose: print(f"GromacsSim::Analyze")

        for sd in self.seeds: sd.Analyze()

        # Do the block averaging
        self.BlockAverage()

        # Categorize the final state of the seed by tilt
        self.CategorizeSeeds()

        print(f"---- Simulation {self.name} analyzed ----")

        if self.verbose: print(f"GromacsSim::Analyze return")

    def WriteData(self):
        r""" Write the data for this simuilation to H5 files
        """
        if self.verbose: print(f"GromacsSim::WriteData")

        # A couple of quick checks to make sure data directories are setup correctly
        if not self.opts.datadir: self.opts.datadir = create_datadir(self.opts.workdir)
        if self.opts.workdir != self.sim_path:
            self.sim_datadir = create_datadir(self.opts.datadir, datadir_name = "{}_data".format(self.name))
        else:
            self.sim_datadir = self.opts.datadir

        #self.final_df.to_hdf(os.path.join(self.sim_datadir, "{}_final.h5".format(self.name)), key = "final_df", mode = "w")
        if self.verbose: print(f"Block average data")
        if self.verbose: print(self.block_df)
        self.block_df.to_hdf(os.path.join(self.sim_datadir, "{}_blocks.h5".format(self.name)), key = "block_df", mode = "w")

        if self.verbose: print(f"GromacsSim::WriteData return")

    def CategorizeSeeds(self):
        r""" Categorize the seeds by tilt
        """
        if self.verbose: print(f"GromacsSim::CategorizeSeeds")

        print(f"WARNING: Categorize Seeds was done in a hurray, check later!")
        print(f"WARNING: Categorize Seeds was done in a hurray, check later!")
        print(f"WARNING: Categorize Seeds was done in a hurray, check later!")

        tilt_criteria = 50.0

        for sd in self.seeds:
            print(f"Categorizing seed: {sd.name}")
           
            target_name = "{}_tilt_mean".format(sd.name)
            tilt_df = self.block_df.filter(regex = target_name)
            final_global_tilt = tilt_df.iloc[-1].to_numpy().flatten()
            is_tilted = (final_global_tilt < tilt_criteria) or (final_global_tilt > 180.0 - tilt_criteria)
            print(f"  Criteria: {is_tilted}, Final tilt: {final_global_tilt}")

        if self.verbose: print(f"GromacsSim::CategorizeSeeds return")

    def BlockAverage(self):
        r""" Do the block averaging of the seeds
        """
        if self.verbose: print(f"GromacsSim::BlockAverage")

        # Dummy graph for us
        fig, ax = plt.subplots(1, 1, figsize = (16,9))

        self.block_dfs = []

        # Loop over seeds
        for sd in self.seeds:
            self.block_average[sd.name] = {}

            # Now do the measurements
            for key,val in self.graph_perseed.items():
                self.block_average[sd.name][key] = {}
                graph = val

                # Timepoints needed now
                timepoints_traj = sd.master_time_df.index/1000.0
                timepoints_forces = sd.master_forces_df.index/1000.0

                if key == "zforce" or key == "perptorque":
                    timepoints = timepoints_forces
                else:
                    timepoints = timepoints_traj

                nblocks = timepoints[-1]/(self.block_size / 1000.0)

                # Get the actual data out from the graph functions
                if self.verbose:
                    print(f"Attempting graph {graph} on sd.label {sd.label}")
                [ylow, yhi, yarr] = graph(sd, ax)
                yarr_np = np.array(yarr)

                # Set up the hi/low parameter choices
                self.ylow_dict[key] = ylow
                self.yhi_dict[key] = yhi

                # Create placeholders for mean and std
                self.block_average[sd.name][key]["mean"] = {}
                self.block_average[sd.name][key]["std"] = {}

                # Loop over the blocks
                for iblock in np.arange(nblocks):
                    start_time = iblock * self.block_size / 1000.0
                    end_time = (iblock + 1) * self.block_size / 1000.0

                    flat_times = timepoints.to_numpy().flatten()
                    time_indices = (flat_times >= start_time) & (flat_times <= end_time)

                    # Get the actual array in the middle of other arrays
                    if yarr_np.ndim == 2:
                        if key == "zforce" or key == "perptorque":
                            mnparr = yarr_np[0,time_indices]
                        else:
                            mnparr = yarr_np[2,time_indices]
                    else:
                        mnparr = yarr_np[time_indices]

                    self.block_average[sd.name][key]["mean"][start_time] = np.mean(mnparr)
                    self.block_average[sd.name][key]["std"][start_time] = np.std(mnparr, ddof=1)

                # Also keep track of this as a dataframe for ease of use
                column_name_mean = "{}_{}_mean".format(sd.name, key)
                column_name_std  = "{}_{}_std".format(sd.name, key)
                index_vals = self.block_average[sd.name][key]["mean"].keys()
                column_vals_mean = self.block_average[sd.name][key]["mean"].values()
                column_vals_std  = self.block_average[sd.name][key]["std"].values()
                column_dict = {column_name_mean: column_vals_mean,
                               column_name_std: column_vals_std}
                df_key = pd.DataFrame(column_dict, index = index_vals)
                self.block_dfs.append(df_key)


        # Clean up the graphs
        plt.close()
        gc.collect()

        # Put the information together
        self.block_df = pd.concat(self.block_dfs, axis=1)

        if self.verbose: print(f"GromacsSim::BlockAverage return")

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

            axarr_mean.set_ylabel(self.named_graphs[key])

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

            #if key == "zforce" or key == "perptorque":
            #    xdata = np.array(self.final_data[key])[:,0]
            #elif key == "zpos":
            #    xdata = np.array(self.final_data[key])[:,2]
            #else:
            #    xdata = np.array(self.final_data[key])

            #plt.figure(fig_dist)
            #fig_range = (self.ylow_dict[key], self.yhi_dict[key])
            #(hist, bin_mids, npoints) = graph_sim_histogram(xdata, ax_dist, mtitle = key, ytitle = "Count", xtitle = self.named_graphs[key], xrange = fig_range)
            #fig_dist.tight_layout()
            #plt.savefig("{}/{}_distribution_{}.pdf".format(self.individual_dir, key, self.name), dpi = fig_dist.dpi)

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

            # Get the information out of the dataframe
            target_xdata_mean   = ".*_{}_mean".format(key1)
            target_xdata_std    = ".*_{}_std".format(key1)
            xdata_mean_df       = self.block_df.filter(regex = target_xdata_mean)
            xdata_std_df        = self.block_df.filter(regex = target_xdata_std)

            target_xdata_mean   = ".*_{}_mean".format(key2)
            target_xdata_std    = ".*_{}_std".format(key2)
            ydata_mean_df       = self.block_df.filter(regex = target_xdata_mean)
            ydata_std_df        = self.block_df.filter(regex = target_xdata_std)

            xdata_mean = xdata_mean_df.iloc[-1,:]
            ydata_mean = ydata_mean_df.iloc[-1,:]
            xdata_std = xdata_std_df.iloc[-1,:]
            ydata_std = ydata_std_df.iloc[-1,:]

            plt.figure(fig_scatter)
            ax_scatter.scatter(x = xdata_mean, y = ydata_mean, zorder = 100, s = 30, marker = 's', color = 'b', facecolors = 'none')
            ax_scatter.errorbar(x = xdata_mean, y = ydata_mean, xerr = xdata_std, yerr = ydata_std, ecolor = 'b', elinewidth = 2, capsize = 5, capthick = 1, zorder = 0, fmt = 'none', marker = 's')

            # If we have stability and the combination is zpos/tilt, then we can graph this
            if key1 == "zpos" and key2 == "tilt" and self.stability_filename:
                df_stability = pd.read_csv(self.stability_filename)
                df_stability.drop_duplicates()

                # Cyan = force untable, magenta = torque unstable, black = stable
                stable_colors = []
                for index,row in df_stability.iterrows():
                    ccolor = 'k'
                    if row['Fstab'] == 0 and row['Tstab'] == 0:
                        ccolor = 'r'
                    elif row['Fstab'] == 0:
                        ccolor = 'c'
                    elif row['Tstab'] == 0:
                        ccolor = 'm'
                    stable_colors.append(ccolor)
                
                predict_scatter = ax_scatter.scatter(df_stability['Z']*10.0, df_stability['Theta'], s = 50, marker = 'o', color = stable_colors, facecolors = stable_colors)
                #predict_legend = ax_scatter.legend(*predict_scatter.legend_elements())

            ax_scatter.set_xlim([self.ylow_dict[key1], self.yhi_dict[key1]])
            ax_scatter.set_ylim([self.ylow_dict[key2], self.yhi_dict[key2]])
            ax_scatter.set_xlabel(self.named_graphs[key1])
            ax_scatter.set_ylabel(self.named_graphs[key2])

            fig_scatter.tight_layout()
            plt.savefig("{}/{}_{}_{}.pdf".format(combodir, key1, key2, self.name), dpi = fig_scatter.dpi)

            # clean up
            plt.close(fig_scatter)
            gc.collect()

        if self.verbose: print(f"GromacsSim::GraphFinalScatter return")
