# XXX: Copyright info here

"""Base class for single simulation object"""

import os
import re
import yaml

import numpy as np

class SimulationBase(object):
    def __init__(self, opts = None):
        print(f"SimulationBase::init")
        self.opts = opts

        # Set up the current working directory properly
        self.cwd = os.getcwd()
        self.ReadOpts()

        self.path = os.path.abspath(self.opts.workdir)
        
        self.yaml_filename = self.opts.yaml
        self.default_yaml = self.GetYamlDict(self.yaml_filename)

        self.ReadData()
        print(f"SimulationBase::init return")

    def ReadOpts(self):
        if not self.opts.workdir:
            self.opts.workdir = os.path.abspath(self.cwd)
        elif not os.path.exists(self.opts.workdir):
            raise IOError("Working directory {} does not exist.".format(
                self.opts.workdir) )
        else:
            self.opts.workdir = os.path.abspath(self.opts.workdir)

    def ReadData(self):
        print(f"SimulationBase::ReadData")
        # Simulation parameters
        self.kT                 = np.float64(self.default_yaml['simulation']['kT'])
        self.bead_size          = np.float64(self.default_yaml['simulation']['bead_size'])
        self.deltatau           = np.float64(self.default_yaml['simulation']['deltatau'])
        self.nsteps             = np.int32(np.float64(self.default_yaml['simulation']['nsteps']))
        self.nwrite             = np.int32(np.float64(self.default_yaml['simulation']['nwrite']))
        self.lbox               = np.float64(self.default_yaml['simulation']['lbox'])
        self.nseed              = np.int32(np.float64(self.default_yaml['simulation']['seed']))
        self.compute_mode       = self.default_yaml['simulation']['mode']
        self.trajectory_file    = self.default_yaml['simulation']['trajectory_file']
        self.init_type          = self.default_yaml['simulation']['init_type']
        self.integrator         = self.default_yaml['simulation']['integrator']

        print(f"SimulationBase::ReadData return")

    def GetYamlDict(self, file_name):
        file_path = os.path.join(self.path, file_name)
        file_dict = ''
        with open(file_path, 'r') as stream: file_dict = yaml.safe_load(stream)
        return file_dict
