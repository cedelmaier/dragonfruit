# XXX: Put a license here

"""Class for a single septin simulation"""

import gc
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from septin_seed import SeptinSeed
from sim_base import SimulationBase

from seed_graph_funcs import *
from sim_graph_funcs import *

class SeptinSim(SimulationBase):
    def __init__(self, path, opts, seedType=SeptinSeed):
        print(f"SeptinSimulation::__init__")
        SimulationBase.__init__(self, path, opts, seedType=seedType)

        # Define functions to graph we want on a per-seed basis
        self.graph_perseed_functions = [ graph_seed_temperature,
                                         graph_seed_pressure,
                                         graph_seed_area ]

        # Now for anything that can be averaged
        self.graph_distribution_functions = [ (graph_sim_membranemodes, 1, 1, (6, 4)) ]

        print(f"SeptinSimulation::__init__ return")

    def Analyze(self):
        r""" Analyze all of the underlying seeds
        """
        for sd in self.seeds: sd.Analyze()
        print(f"- Sim {self.name} analyzed -")

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
        self.GraphDistributions()

    def GraphDynamicData(self):
        r""" Graph dynamic data (similar to seed methods) for simulation
        """
        fig, axarr = plt.subplots(3, 2, figsize=(15,10))
        colors = mpl.cm.rainbow(np.linspace(0,1,len(self.seeds)))

        # Zip together the list like we had before to be clever
        # This first loop zips together the functions we want to graph and the axarr we generated
        for graph, axr in zip(self.graph_perseed_functions, axarr):
            # The second loop creates different colors for each seed
            for sd, col in zip(self.seeds, colors):
                yarr = graph(sd, axr[0], color=col, xlabel = False) # Actual magic of the graphing call

        #axarr[1,0].legend(loc = 'center left', bbox_to_anchor=(2.2,-0.19))
        fig.tight_layout()
        plt.savefig("{}_{}.pdf".format(os.path.join(self.sim_datadir, 'dynamicdata'), self.name), dpi=fig.dpi)
        
        # Clean up
        fig.clf()
        plt.close()
        gc.collect()
        mpl.rcdefaults()

    def GraphDistributions(self):
        r""" Graph distribution data for simulation
        """
        # This is a 'clever' loop over the correct distributions as defined earlier

        for gf, rn, cn, fs in self.graph_distribution_functions:
            fig, axarr = plt.subplots(rn, cn, figsize = fs)
            gf(self, axarr)
            plt.figure(fig.number)
            fig.tight_layout()
            pruned_name = self.name[0:127]
            fig.savefig("{}_{}.pdf".format(os.path.join(self.sim_datadir, gf.__name__), pruned_name), dpi=fig.dpi)

            # Clean up
            fig.clf()
            plt.close()
            gc.collect()

        mpl.rcdefaults()

