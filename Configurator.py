# Configuration class for Lipid systems (extensible in future)
import argparse
import datetime
import itertools
import os
import sys
import yaml

import numpy as np

# Class definition
class Configurator(object):
    def __init__(self, opts):
        self.default_yaml = self.GetYamlDict(opts.default_file)

    # Get a YAML dictionary for the file
    def GetYamlDict(self, filename):
        file_dict = ''
        with open(filename, 'r') as stream: file_dict = yaml.safe_load(stream)
        return file_dict

