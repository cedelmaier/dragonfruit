# Configuration class for Lipid systems (extensible in future)
import hoomd
import hoomd.md as md
import gsd.hoomd

import itertools
import os
import random
import sys
import yaml

import numpy as np

from common import *
#from ahdomain import ahdomain
from block_ahdomain import block_ahdomain
from helix_ahdomain import helix_ahdomain
from membrane import membrane

# Configurator class definition
class Configurator(object):
    def __init__(self, opts):
        self.default_yaml = self.GetYamlDict(opts.default_file)

        self.SetDefaults()

    def SetDefaults(self):
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

        # Set the internal data classes
        self.lipids = None
        self.ahdomain = None

        # Set the random seed
        random.seed(self.nseed)

        # Check integrator settings to make sure it's okay!
        if (self.integrator != 'langevin' and self.integrator != 'NPT' and self.integrator != 'NPH'):
            print(f"Integrator {self.integrator} not used for this system, exiting!")
            sys.exit(1)

        # Unfortunately, now have some gymnastics for checking if something exists and setting
        if 'pdamp' in self.default_yaml['simulation']:
            self.pdamp = np.float64(self.default_yaml['simulation']['pdamp'])
        else:
            self.pdamp = None

        if 'tau' in self.default_yaml['simulation']:
            self.tau = np.float64(self.default_yaml['simulation']['tau'])
        else:
            self.tau = None

        if 'tauS' in self.default_yaml['simulation']:
            self.tauS = np.float64(self.default_yaml['simulation']['tauS'])
        else:
            self.tauS = None

        # Are we reading in previous information?
        if self.init_type == 'all':
            self.init_filename = ''
        elif self.init_type == 'read_gsd':
            self.init_filename = self.default_yaml['simulation']['init_filename']
        else:
            print(f"Need to specify a correct initialization type, tried {self.init_type}, exiting!")
            sys.exit(1)

        # Lipid parameters
        if 'membrane' in self.default_yaml:
            self.lipids = membrane(self.bead_size, self.default_yaml)

        # AH parameters
        if 'ah_domain' in self.default_yaml:
            self.ahtype = self.default_yaml['ah_domain']['polymer_type']
            if self.ahtype == 'block_copolymer':
                self.ahdomain = block_ahdomain(self.bead_size, self.default_yaml)
            elif self.ahtype == 'helix_block_copolymer':
                self.ahdomain = helix_ahdomain(self.bead_size, self.default_yaml)
            else:
                print(f"AH domain type {self.ahtype} not implemented, exiting!")
                sys.exit(1)

    # Print information
    def PrintInformation(self, snap):
        print(f"--------")
        print(f"System information")
        print(f"Compute mode            = {self.compute_mode}")
        print(f"Integrator              = {self.integrator}")
        print(f"Delta tau               = {self.deltatau}")
        print(f"Simulation time (tau)   = {self.deltatau*self.nsteps}")
        print(f"kBT                     = {self.kT}")
        print(f"seed                    = {self.nseed}")
        print(f"box                     = ({snap.configuration.box[0]}, {snap.configuration.box[1]}, {snap.configuration.box[2]})")
        # Configurable parameters
        if self.tau: print(f"tau                     = {self.tau}")
        if self.pdamp: print(f"pdamp                   = {self.pdamp}")
        if self.tauS: print(f"tauS                    = {self.tauS}")

        if self.lipids:
            self.lipids.PrintInformation(snap)
        if self.ahdomain:
            self.ahdomain.PrintInformation(snap)

    # Initialize information
    def Init(self, snap):
        # Create a default configuration box size if we don't have a membrane
        snap.configuration.box = [self.lbox, self.lbox, self.lbox, 0, 0, 0]

        if self.lipids:
            self.lipids.InitMembrane(snap)
        if self.ahdomain:
            self.ahdomain.InitAH(snap)

    # Get a YAML dictionary for the file
    def GetYamlDict(self, filename):
        file_dict = ''
        with open(filename, 'r') as stream: file_dict = yaml.safe_load(stream)
        return file_dict

