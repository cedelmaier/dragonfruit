# XXX: Copyright info here

"""Base class for simulation seed analysis"""

import os
import re
import yaml

class SeedBase(object):
    def __init__(self, path, opts = None):
        self.opts = opts
        self.path = os.path.abspath(path)
        self.name = self.path.split('/')[-1]
        self.label, self.seed_num = self.MakeSeedLabel()
        
        self.yaml_filename = opts.yaml
        self.default_yaml = self.GetYamlDict(self.yaml_filename)

    def MakeSeedLabel(self):
        snum = re.findall("\d*\.?\d+", self.name)
        try:
            slabel = r's $=$ {}'.format(snum[-1])
            return slabel, int(snum[-1])
        except:
            print("Could not make seed label. Using name of seed directory.")
            return self.name, 0

    def GetYamlDict(self, file_name):
        file_path = os.path.join(self.path, file_name)
        file_dict = ''
        with open(file_path, 'r') as stream: file_dict = yaml.safe_load(stream)
        return file_dict


