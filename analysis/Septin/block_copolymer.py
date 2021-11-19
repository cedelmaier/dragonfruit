# XXX: Put a license here

"""Helper class for AH domain block copolymers"""

import numpy as np
import pandas as pd

class BlockCopolymer(object):
    def __init__(self, yaml_dict, ah_start_nx):
        r""" Initialize a block copolymer helper

        Note, the ah_start_nx tells us the starting particle index of the AH domain
        """
        self.default_yaml = yaml_dict
        self.ah_start_nx = ah_start_nx
        
        self.nbeads = self.default_yaml['ah_domain']['nbeads']
        self.nrepeat = self.default_yaml['ah_domain']['nrepeat']
