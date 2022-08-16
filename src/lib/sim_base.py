# XXX: Copyright info here

"""Base class for single simulation object"""

import ast
import os
import re
import yaml

import numpy as np

from seed_base import SeedBase

class SimulationBase(object):
    def __init__(self, sim_path, opts, seedType = SeedBase):
        print(f"SimulationBase::__init__")
        
        self.sim_path = os.path.abspath(sim_path)
        self.name = self.sim_path.split('/')[-1]
        self.title = self.MakeSimTitle(self, sim_path)
        self.params = self.MakeSimParamDict(self, sim_path)
        self.opts = opts

        self.seedType = seedType
        self.seeds = []

        self.CollectSeeds(seedType)

        print(f"SimulationBase::__init__ return")

    @staticmethod
    def MakeSimTitle(self, sim_path, uc = ''):
        r""" Create a simulation title based on the parameters
        """
        name = sim_path.split('/')[-1]
        # Check for an MD5 hash to make sure we don't bomb
        is_hash = re.findall(r"([a-fA-F\d]{32})", name) 
        if len(is_hash) > 0:
            label = r"hash {}".format(is_hash)
            return label
        param = re.findall("[a-zA-Z]+", name)
        p_value = re.findall("\d*\.?\d+", name)

        # Specific check for scientific notation
        while 'e' in param:
            param.remove('e')

        match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
        p_value2 = re.findall(match_number, name)

        label = r''
        slabel = r'{} $=$ {}{}, '
        for p, v in zip(param, p_value2):
            if uc:
                label += slabel.format(p,
                        str(ast.literal_eval(v)*uc[p][1]), uc[p][0])
            else:
                label += slabel.format(p,
                        str(ast.literal_eval(v)), '')
        return label[:-2]

    @staticmethod
    def MakeSimParamDict(self, sim_path, uc=''):
        r""" Create the parameter list that this simulation uses
        """
        params = {}
        sim_name = sim_path.split('/')[-1]
        # Check for an MD5 hash to make sure we don't bomb
        is_hash = re.findall(r"([a-fA-F\d]{32})", sim_name)
        if len(is_hash) > 0:
            return params
        p_name = re.findall("[a-zA-Z]+", sim_name)
        p_value = re.findall("\d*\.?\d+", sim_name)

        # Check for hash values as well
        while 'e' in p_name:
            p_name.remove('e')
        match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
        p_value2 = re.findall(match_number, sim_name)

        for p, v in zip(p_name, p_value2):
            if uc:
                params[p] = ast.literal_eval(v)*uc[p][1]
            else:
                params[p] = ast.literal_eval(v)

        return params

    def CollectSeeds(self, seedType):
        r""" Collect the underlying seeds for this simulation
        """
        # Add a check for the gromacs type numbering system
        seed_pattern = re.compile("[s|N]\d+")

        # Make list of seeds of type seedType
        self.seeds = [seedType(os.path.join(self.sim_path,sd), self.opts)
                for sd in next(os.walk(self.sim_path))[1] if seed_pattern.match(sd)]
        self.seeds.sort(key=lambda sd: sd.seed_num)

    def GraphSimulation(self):
        raise NotImplementedError()
