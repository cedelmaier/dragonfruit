# Configuation class for ah domains in HOOMD
# XXX: Maybe make a dataclass in the future
import hoomd
import hoomd.md as md
import gsd.hoomd

import itertools
import os
import random
import sys
import yaml

import numpy as np

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from common import *

class ahdomain(object):
    def __init__(self, bead_size, yaml_file):
        # AH-domain parameters
        self.nah              = np.int32(np.float64(yaml_file['ah_domain']['n_ah']))
        self.polymer_type     = yaml_file['ah_domain']['polymer_type']
        self.gamma            = np.float64(yaml_file['ah_domain']['gamma'])
        self.mass             = np.float64(yaml_file['ah_domain']['mass'])
        self.A                = np.float64(yaml_file['ah_domain']['A'])
        self.Bsurface         = np.float64(yaml_file['ah_domain']['B_surface'])
        self.Bintermediate    = np.float64(yaml_file['ah_domain']['B_intermediate'])
        self.Bdeep            = np.float64(yaml_file['ah_domain']['B_deep'])
        #self.Bself            = np.float64(yaml_file['ah_domain']['B_self'])
        self.Bself_hydrophobic  = np.float64(yaml_file['ah_domain']['B_self_hydrophobic'])
        self.Bself_other        = np.float64(yaml_file['ah_domain']['B_self_other'])
        self.kbond            = np.float64(yaml_file['ah_domain']['kbond'])
        self.kbend            = np.float64(yaml_file['ah_domain']['kbend'])
        self.is_init          = yaml_file['ah_domain']['is_init']

    # Initialize the AH domains
    def InitAH(self, snap):
        if self.is_init:
            print(f"Already think AH is initialized")
        else:
            self.CreateAH(snap)
        self.is_init = True

    # Create some number of block copolyerms in the simulation
    def CreateAH(self, snap):
        raise NotImplementedError("Please Implement this method")
       
    def PrintInformation(self, snap):
        raise NotImplementedError("Please Implement this method")
