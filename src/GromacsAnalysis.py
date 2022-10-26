#!/usr/bin/env python3

# XXX: Put a license here

"""Main analysis script for membranes with AH domains"""

import argparse
import logging
import pickle
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
sys.path.append(os.path.join(os.path.dirname(__file__), 'gromacs'))
from common import create_datadir
from stylelib.common_styles import *
from gromacs_seed import GromacsSeed
from gromacs_sim import GromacsSim
from gromacs_run import GromacsRun

def parse_args():
    parser = argparse.ArgumentParser(prog='GromacsAnalysis.py')

    # General options
    parser.add_argument('-sd', '--seed', action = 'store_true',
            help = 'Run analysis on a single seed')
    parser.add_argument('-sim', '--simulation', action = 'store_true',
            help = 'Run analysis on a single simulation (collection of seeds)')

    parser.add_argument('-d', '--workdir', action = 'store_true',
            help = 'Working directory')
    parser.add_argument('--datadir', type=str,
            help = 'Data directory where analyzed files will be written to. Defaults to {workdir}/data')
    parser.add_argument('--yaml', type = str,
            help = 'YAML file to read')

    # Flow control options
    parser.add_argument('-A', '--analysis', action = 'store_true',
            help = 'Analyze data from MD simulation')
    parser.add_argument('-G', '--graph', action = 'store_true',
            help = 'Graph data from MD simulation')
    parser.add_argument('-W', '--write', action = 'store_true',
            help = 'Write data to HD5 and pickle files')
    parser.add_argument('-F', '--force', action = 'store_true',
            help = 'Force complete analysis of simulation(s)')

    # Add verbosity control
    parser.add_argument('-v', '--verbose', action="store_true",
            help = 'Verbose output')
    parser.add_argument('-t', '--trace', action ="store_true",
            help = 'Trace output')

    opts = parser.parse_args()
    return opts

class GromacsAnalysis(object):
    r""" Gromacs analysis
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
            raise IOError("Working directory {} does not exist.".format(self.opts.workdir))
        else:
            self.opts.workdir = os.path.abspath(self.opts.workdir)

        self.verbose = False
        if self.opts.verbose:
            self.verbose = True

        # If there is some central amount of data, create the directory
        if (self.opts.write or self.opts.graph):
            self.opts.datadir = create_datadir(self.opts.workdir)

    def ProgOpts(self):
        r""" Run selected commands
        """
        
        if self.opts.seed:
            self.AnalyzeSeed()
        elif self.opts.simulation:
            self.AnalyzeSimulation()
        else:
            self.AnalyzeRun()

    def AnalyzeSeed(self):
        r""" Analyze a single simulation seed
        """
        sd = GromacsSeed(self.opts.workdir, self.opts)
        sd.Analyze()

        if self.opts.graph:
            plt.style.use(septin_runs_stl)

            # Plot the Z distance versus time
            fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
            sd.GraphZdist(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_zdist.pdf', dpi=fig.dpi)

            # Try doing a plot with the lipid head groups included
            fig, axarr = plt.subplots(1, 1, figsize=(15,10))
            sd.GraphZpos(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_zpos.pdf', dpi=fig.dpi)

            # Plot the helicity of the protein
            fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
            sd.GraphHelicity(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_helicity.pdf', dpi=fig.dpi)

            ## Plot helix information
            #fig, axarr = plt.subplots(2, 1, figsize = (15, 10))
            #sd.graph_helix_analysis(axarr)
            #fig.tight_layout()
            #fig.savefig('gromacs_helix_analysis.pdf', dpi=fig.dpi)

            # Plot the global tilts
            fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
            sd.GraphGlobalTilt(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_global_tilt.pdf', dpi=fig.dpi)

            # Plot the mean location of the upper leaflet
            fig, axarr = plt.subplots(1, 2, figsize = (15, 10))
            sd.graph_avg_z_surface(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_avg_zsurf.pdf', dpi=fig.dpi)

            # Plot the average mean curvature
            fig, axarr = plt.subplots(1, 2, figsize = (15, 10))
            sd.graph_avg_mean_curvature(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_avg_meancurv.pdf', dpi = fig.dpi)

            # Plot the average gaussian curvature
            fig, axarr = plt.subplots(1, 2, figsize = (15, 10))
            sd.graph_avg_gaussian_curvature(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_avg_gausscurv.pdf', dpi = fig.dpi)

            # Plot the angle of the dipole moment with the Z-axis
            fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
            sd.GraphPDipoleTilt(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_pdipole_tilt.pdf', dpi=fig.dpi)

            # Plot the angle of the dipole moment with the Z-axis
            fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
            sd.GraphHelixPDipoleAngle(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_helixpdipole_angle.pdf', dpi=fig.dpi)

            # Plot the magnitude of the electric dipole moment
            fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
            sd.GraphPDipoleMoment(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_pdipolemoment.pdf', dpi = fig.dpi)

            # Plot the Z component of the force
            fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
            sd.GraphZForce(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_zforce.pdf', dpi = fig.dpi)

            # Plot the perpendicular torque to the helix
            fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
            sd.GraphPerpTorque(axarr)
            fig.tight_layout()
            fig.savefig('gromacs_perptorque.pdf', dpi = fig.dpi)

    def AnalyzeSimulation(self):
        r""" Analyze a simulation (collection of seeds)
        """
        sim = GromacsSim(self.opts.workdir, opts)
        sim.Analyze()
        if self.opts.graph:
            plt.style.use(septin_poster_stl)
            sim.GraphSimulation()

        if self.opts.write:
            sim.WriteData()

    def AnalyzeRun(self):
        r""" Analyze a run (collection of simulations)
        """
        run = GromacsRun(opts, self.opts.workdir)

        if self.opts.graph:
            plt.style.use(septin_poster_stl)
            run.GraphRun()

##########################################
if __name__ == "__main__":
    opts = parse_args()
    x = GromacsAnalysis(opts)
