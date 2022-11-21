# General submission script for mucus simulations
import hoomd
import hoomd.md as md
import gsd.hoomd

import argparse
import datetime
import itertools
import os
import sys
import yaml

import numpy as np

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'mucus'))
from mucus_seed import MucusSeed

# Create a rudimentary parser to get the number of steps, the write frequency,
# and the box length
def parse_args():
    parser = argparse.ArgumentParser(prog='MucusVolta.py')

    # General options
    parser.add_argument('--yaml', type=str, default='mucus.default.yaml',
            help='YAML configuration file')

    parser.add_argument('-d', '--workdir', action = 'store_true',
            help = 'Working directory')

    # Add verbosity control
    parser.add_argument('-v', '--verbose', action="store_true",
                        help = 'Verbose output')

    opts = parser.parse_args()

    return opts

# Create a status line maker for our output
class Status():

    def __init__(self, sim):
        self.sim = sim

    @property
    def seconds_remaining(self):
        try:
            return (self.sim.final_timestep - self.sim.timestep) / self.sim.tps
        except ZeroDivisionError:
            return 0

    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining))

###############################################################################
# Main program start
###############################################################################
if __name__ == "__main__":
    # Parse those args
    opts = parse_args()

    cwd = os.getcwd()
    if not opts.workdir:
        opts.workdir = os.path.abspath(cwd)
    elif not os.path.exists(opts.workdir):
        raise IOError("Working directory {} does not exist, exiting.".format(
            opts.workdir))
    else:
        opts.workdir = os.path.abspath(opts.workdir)

    # Configurator is just a septin seed, since we moved the functionality into that class
    configurator = MucusSeed(opts.workdir, opts)

    configurator.Configure()

    configurator.PrintInformation()
