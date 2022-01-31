#!/usr/bin/env python3

# XXX: Put a license here

"""Main analysis scripts for septins and membranes."""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get the proper library directories
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
sys.path.append(os.path.join(os.path.dirname(__file__), 'septins'))
from stylelib.common_styles import septin_runs_stl
from septin_seed import SeptinSeed

def parse_args():
    parser = argparse.ArgumentParser(prog='SeptinAnalysis.py')

    # General options
    parser.add_argument('-sd', '--seed', action = 'store_true',
            help = 'Run analysis on a single seed')

    parser.add_argument('-d', '--workdir', action = 'store_true',
            help = 'Working directory')

    parser.add_argument('--yaml', type = str,
            help = 'YAML file to read')

    # Control options
    parser.add_argument('-A', '--analyze', action = 'store_true',
            help = 'Analyze data from simulation(s)')

    parser.add_argument('-F', '--force', action = 'store_true',
            help = 'Force complete analysis of simulation(s)')

    parser.add_argument('-G', '--graph', action = 'store_true',
            help = 'Graph data after analysis has been performed.')

    opts = parser.parse_args()
    return opts

class SeptinAnalysis(object):
    r"""Septin and membrane analysis
    """
    def __init__(self, opts):
        self.opts = opts
        self.cwd = os.getcwd()

        self.ReadOpts()

        self.ProgOpts()

    def ReadOpts(self):
        if not self.opts.workdir:
            self.opts.workdir = os.path.abspath(self.cwd)
        elif not os.path.exists(self.opts.workdir):
            raise IOError("Working directory {} does not exist.".format(
                self.opts.workdir) )
        else:
            self.opts.workdir = os.path.abspath(self.opts.workdir)

    def ProgOpts(self):
        r"""Run selected commands
        """

        if self.opts.seed:
            self.AnalyzeSeed()

    def AnalyzeSeed(self):
        r"""Analyze a single simulation seed
        """
        sd = SeptinSeed(self.opts.workdir, self.opts)
        sd.Analyze()

        if self.opts.graph:
            plt.style.use(septin_runs_stl)

            # Plot dynamic variables
            fig, axarr = plt.subplots(3, 1, figsize=(15, 10))
            sd.GraphDynamic(axarr)
            fig.tight_layout()
            fig.savefig('septin_parameters.pdf', dpi=fig.dpi)

            # Plot diffusion if it exists
            fig, axarr = plt.subplots(1, 1, figsize=(15,10))
            sd.GraphSimpleMSD(axarr)
            fig.tight_layout()
            fig.savefig('septin_simple_msd.pdf', dpi=fig.dpi)

            # Plot distribution data
            fig, axarr = plt.subplots(3, 2, figsize=(15,10))
            sd.GraphDistributions(axarr)
            fig.tight_layout()
            fig.savefig('septin_distributions.pdf', dpi=fig.dpi)

            # Plot membrane modes
            fig, axarr = plt.subplots(1, 1, figsize=(15,10))
            sd.GraphMembraneModes(axarr)
            fig.tight_layout()
            fig.savefig('septin_membranemodes.pdf', dpi=fig.dpi)
        


##########################################
if __name__ == "__main__":
    opts = parse_args()
    x = SeptinAnalysis(opts)
