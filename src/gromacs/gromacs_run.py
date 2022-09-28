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
        self.main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data')

        # Wrap up stuff into simulation packages
        self.simulation_packages = {
                               'polarmonomer_aglipid_rotxN_50mMKCl':  [[
                                                                      'unbiased_zdepth00_rotx0_helix_50mMKCl',
                                                                      'unbiased_zdepth00_rotx90_helix_50mMKCl',
                                                                      'unbiased_zdepth00_rotx180_helix_50mMKCl',
                                                                      'unbiased_zdepth00_rotx270_helix_50mMKCl',
                                                                      ],
                                                                      [
                                                                      r'AGMAGL 0 deg 50 mM KCl',
                                                                      r'AGMAGL 90 deg 50 mM KCl',
                                                                      r'AGMAGL 180 deg 50 mM KCl',
                                                                      r'AGMAGL 270 deg 50 mM KCl',
                                                                      ]],
                               'polarmonomer_aglipid_rotxNtrans_50mMKCl':  [[
                                                                      'unbiased_zdepth00_rotx0_helix_50mMKCl',
                                                                      'unbiased_zdepth00_rotx90_helix_50mMKCl',
                                                                      'unbiased_zdepth00_rotx180_helix_50mMKCl',
                                                                      'unbiased_zdepth00_rotx270_helix_50mMKCl',
                                                                      'unbiased_zdepth00_trans_helix_50mMKCl',
                                                                      ],
                                                                      [
                                                                      r'AGMAGL 0 deg 50 mM KCl',
                                                                      r'AGMAGL 90 deg 50 mM KCl',
                                                                      r'AGMAGL 180 deg 50 mM KCl',
                                                                      r'AGMAGL 270 deg 50 mM KCl',
                                                                      r'AGMAGL Trans 50 mM KCl',
                                                                      ]],
                               'polarmonomer_aglipid_rotxN_150mMKCl':  [[
                                                                      'unbiased_150mMKCl',
                                                                      'agmonomer_11x11_zdepth00_rotx90_150KCl',
                                                                      'agmonomer_aglipid_11x11_zdepth00_rotx180_150mMKCl',
                                                                      'agmonomer_aglipid_11x11_zdepth00_rotx270_150mMKCl',
                                                                      ],
                                                                      [
                                                                      r'AGMAGL 0 deg 150 mM KCl',
                                                                      r'AGMAGL 90 deg 150 mM KCl',
                                                                      r'AGMAGL 180 deg 150 mM KCl',
                                                                      r'AGMAGL 270 deg 150 mM KCl',
                                                                      ]],
                               'neutralmonomer_aglipid_rotxN_50mMKCl':  [[
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx0_50mMKCl',
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx90_50mMKCl',
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx175_50mMKCl',
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx270_50mMKCl',
                                                                      ],
                                                                      [
                                                                      r'RFMAGL 0 deg 50 mM KCl',
                                                                      r'RFMAGL 90 deg 50 mM KCl',
                                                                      r'RFMAGL 175 deg 50 mM KCl',
                                                                      r'RFMAGL 270 deg 50 mM KCl',
                                                                      ]],
                               'neutralmonomer_aglipid_rotxNall_50mMKCl':  [[
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx0_50mMKCl',
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx90_50mMKCl',
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx175_50mMKCl',
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx265_50mMKCl',
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx270_50mMKCl',
                                                                      ],
                                                                      [
                                                                      r'RFMAGL 0 deg 50 mM KCl',
                                                                      r'RFMAGL 90 deg 50 mM KCl',
                                                                      r'RFMAGL 175 deg 50 mM KCl',
                                                                      r'RFMAGL 265 deg 50 mM KCl',
                                                                      r'RFMAGL 270 deg 50 mM KCl',
                                                                      ]],
                                'neutralmonomer_aglipid_rotxN_150mMKCl': [[
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx0_150KCl',
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx90_150mMKCl',
                                                                      'rfmonomer_aglipid_11x11_zdepth00_rotx270_150mMKCl',
                                                                      ],
                                                                      [
                                                                      r'RFMAGL 0 deg 150 mM KCl',
                                                                      r'RFMAGL 90 deg 150 mM KCl',
                                                                      r'RFMAGL 270 deg 150 mM KCl',
                                                                      ]],
                                                                      }

        # Get the simulations to run for this setup
        self.simulations_to_run = [
                                    'polarmonomer_aglipid_rotxN_50mMKCl',
                                    #'polarmonomer_aglipid_rotxNtrans_50mMKCl',
                                    'polarmonomer_aglipid_rotxN_150mMKCl',
                                    'neutralmonomer_aglipid_rotxN_50mMKCl',
                                    #'neutralmonomer_aglipid_rotxNall_50mMKCl',
                                    'neutralmonomer_aglipid_rotxN_150mMKCl',
                                   ]

        self.simulation_text_data = {}
        self.simulation_text_data['polarmonomer_aglipid_rotxN_50mMKCl'] = 'Polar 50 mM KCl'
        self.simulation_text_data['polarmonomer_aglipid_rotxNtrans_50mMKCl'] = 'Polar 50 mM KCl (with trans)'
        self.simulation_text_data['polarmonomer_aglipid_rotxN_150mMKCl'] = 'Polar 150 mM KCl'
        self.simulation_text_data['neutralmonomer_aglipid_rotxN_50mMKCl'] = 'Neutral 50 mM KCl'
        self.simulation_text_data['neutralmonomer_aglipid_rotxNall_50mMKCl'] = 'Neutral 50 mM KCl (all)'
        self.simulation_text_data['neutralmonomer_aglipid_rotxN_150mMKCl'] = 'Neutral 150 mM KCl'

        # Some final resting places for data, last 25 ns of the simulations
        self.measurement_names = [
                             'zpos',
                             'helix',
                             'tilt',
                             'pdip',
                             'pmom',
                             'hp',
                            ]
        self.measurement_text = {}
        self.measurement_text['zpos'] = (r"Z ($\AA$)", -30.0, 30.0)
        self.measurement_text['helix'] = (r"Helicity (AU)", 0.0, 1.05)
        self.measurement_text['tilt'] = (r"Helix tilt (deg)", 0.0, 180.0)
        self.measurement_text['pdip'] = (r"Helix dipole tilt (deg)", 0.0, 180.0)
        self.measurement_text['pmom'] = (r"Helix electric dipole moment (?)", 0.0, 50.0)
        self.measurement_text['hp'] = (r"Helix-dipole angle (deg)", 0.0, 180.0)

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

        # Common parameter choices for this analysis, report these to the user
        self.common_start_interval  =  50000.0 # First 50 ns
        self.common_end_interval    = 150000.0 # Last 50 ns

        print(f"Beginning of simulation considered to be 0.0 to {self.common_start_interval}")
        print(f"End of simulation considered to be {self.common_end_interval} to end")
        
        if opts.verbose: print(f"GromacsRun::__init__ return")

    def AggregateSimilar(self):
        r""" Aggregate similar condition simulations
        """
        if self.verbose: print(f"GromacsRun::AggregateSimilar")
        # Just need to get the similar simulations from the list in the map, as they are already
        # packaged together. Unfortunately, this does duplicate some of what is in GromacsSim
        # as we are essentially pretending that we have a sim information, even though it comes
        # from multiple places
        print(f"Aggregating similar simulations and graphing")
        for simtorun in self.simulations_to_run:
            print(f"  Simulation block {simtorun}")

            # Try to create a subdirectory if it needs to exist
            simdir = create_datadir(self.opts.datadir, simtorun)

            # How many seeds (total) do we have? Find the simulation, subsimulation, and all the seeds in it
            self.simulations = []
            self.seeds = []
            for subsim in self.simulation_packages[simtorun][0]:
                print(subsim)
                # Load this simulation
                simpath = os.path.join(self.main_path, subsim)
                current_simulation = GromacsSim(simpath, self.opts)
                current_simulation.Analyze()

                self.simulations.append(current_simulation)
                for sd in current_simulation.seeds:
                    self.seeds.append(sd)

            nseeds = len(self.seeds)
            # Create a plot for each variable
            for key,val in self.graph_perseed.items():
                print(f"    Graph: {key}")

                # Create different plots for different slices of the variables
                #fig_start, axarr_start  = plt.subplots(1, 2, figsize = (25,16))
                #fig_end,   axarr_end    = plt.subplots(1, 2, figsize = (25,16))
                fig,       axarr        = plt.subplots(1, 2, figsize = (25,16))

                graph = val #alias the graph
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
                yarr_avg = None
                num_seeds = 0
                for sd, col in zip(self.seeds, colors):
                    [ylow, yhi, yarr] = graph(sd, axarr[0], color=col, xlabel = False, dolegend = dolegend)
                    if yarr_avg is None:
                        yarr_avg = np.zeros_like(yarr)
                    yarr_avg = np.add(yarr_avg, yarr)
                    num_seeds += 1

                # Actually do the average
                yarr_avg /= np.float32(num_seeds)

                # Check for multiple plots
                if yarr_avg.ndim == 2:
                    # Check if we are doing forces/torques, or just the membrane
                    if key == "zforce" or key == "perptorque":
                        axarr[1].plot(timepoints, yarr_avg[0][:], linestyle = "solid", label = "All")
                        axarr[1].plot(timepoints, yarr_avg[1][:], linestyle = "dotted", label = "Electrostatic")
                        axarr[1].plot(timepoints, yarr_avg[2][:], linestyle = "dashed", label = "Non-electrostatic")
                        axarr[1].legend()
                    else:
                        axarr[1].plot(timepoints, yarr_avg[0][:], color = "slategrey")
                        axarr[1].plot(timepoints, yarr_avg[1][:], color = "slategrey")
                        axarr[1].plot(timepoints, yarr_avg[2][:])
                else:
                    axarr[1].plot(timepoints, yarr_avg)
                axarr[1].set_ylim([ylow, yhi])

                # Set axis labels
                axarr[0].set_xlabel(r'Time (ns)')
                axarr[1].set_xlabel(r'Time (ns)')

                fig.tight_layout()
                plt.savefig("{}/{}_{}.pdf".format(simdir, key, simtorun), dpi = fig.dpi)

                # Dirty trick to get the early/late times?
                axarr[0].set_xlim([0.0, self.common_start_interval/1000.0])
                axarr[1].set_xlim([0.0, self.common_start_interval/1000.0])
                fig.tight_layout()
                plt.savefig("{}/{}_{}_early.pdf".format(simdir, key, simtorun), dpi = fig.dpi)

                axarr[0].set_xlim([self.common_end_interval/1000.0, 200000.0/1000.0])
                axarr[1].set_xlim([self.common_end_interval/1000.0, 200000.0/1000.0])
                fig.tight_layout()
                plt.savefig("{}/{}_{}_late.pdf".format(simdir, key, simtorun), dpi = fig.dpi)


                # Clean up
                fig.clf()
                plt.close()
                gc.collect()
                #mpl.rcdefaults()

        if self.verbose: print(f"GromacsRun::AggregateSimilar return")
