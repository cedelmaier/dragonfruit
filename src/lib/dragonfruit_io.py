# Common I/O routines for various specialized file formats

import os
from pathlib import Path
import re
import shlex
import sys

import numpy as np
import pandas as pd

# Read an XVG file (generic) and return it as a dataframe
def read_xvg(fname, identifier = ''):
    r""" Read XVG file legends and data
    
    Returns metadata, raw information, and a pandas dataframe containing all of the information
    """

    _ignored = set(('legend', 'view'))
    _re_series = re.compile('s[0-9]+$')
    _re_xyaxis = re.compile('[xy]axis$')
    
    metadata = {}
    num_data = []
    
    metadata['labels'] = {}
    metadata['labels']['series'] = []
    
    ff_path = os.path.abspath(fname)
    if not os.path.isfile(ff_path):
        raise IOError('File not readable: {0}'.format(ff_path))
    
    with open(ff_path, 'r') as fhandle:
        for line in fhandle:
            line = line.strip()
            if line.startswith('@'):
                tokens = shlex.split(line[1:])
                if tokens[0] in _ignored:
                    continue
                elif tokens[0] == 'TYPE':
                    if tokens[1] != 'xy':
                        raise ValueError('Chart type unsupported: \'{0}\'. Must be \'xy\''.format(tokens[1]))
                elif _re_series.match(tokens[0]):
                    metadata['labels']['series'].append(tokens[-1])
                elif _re_xyaxis.match(tokens[0]):
                    metadata['labels'][tokens[0]] = tokens[-1]
                elif len(tokens) == 2:
                    metadata[tokens[0]] = tokens[1]
                else:
                    print('Unsupported entry: {0} - ignoring'.format(tokens[0]), file=sys.stderr)
            elif line[0].isdigit():
                num_data.append(map(float, line.split()))

    num_data = zip(*num_data)
  
    if not metadata['labels']['series']:
        for series in range(len(num_data) - 1):
            metadata['labels']['series'].append('')
    
    # Hack to get from python2 to python3
    num_data = list(num_data)

    # Create a dataframe to return as well, for future ease, do the fast version
    dfs = []
    for i in range(len(num_data[1:])):
        column_name = metadata['labels']['series'][i].replace(" ","") + "_" + identifier
        dfs.append(pd.DataFrame(num_data[i+1], columns = [column_name], index = num_data[0]))
    df = pd.concat(dfs, axis = 1)
    df.index.name = 'Time(ps)'

    return metadata, num_data, df

# ITP Parser class
class dfruit_ITPParser(object):
    def __init__(self, path, verbose = False):
        r""" Create a new ITP Parser object
        """
        self.path = path
        self.verbose = verbose
        if self.verbose: print(f"dfruit_ITPParser::__init__")
        print(f"{path}")
        if self.verbose: print(f"dfruit_ITPParser::__init__ return")

    def parse(self):
        r""" Parser the ITP topology

        Based off of the MDAnalysis ITP parser, but customized for our needs
        https://docs.mdanalysis.org/2.1.0/_modules/MDAnalysis/topology/ITPParser.html#ITPParser

        Returns
        --------
        Loaded topology information for this file
        """
        print(f"ERROR: Experimental feature, exiting!")
        sys.exit(1)
        if self.verbose: print(f"dfruit_ITPParser::parse")

        # Initial parser
        self.parser = self._pass 

        with open(self.path, 'r') as itpfile:
            for line in itpfile:
                print(line.rstrip())
                if '[' in line and ']' in line:
                    section = line.split('[')[1].split(']')[0].strip()

                    if section == 'moleculetype':
                        self.parser = self.parser_moleculetype
                    else:
                        self.parser = self._pass

                else:
                    self.parser(line)







        if self.verbose: print(f"dfruit_ITPParser::parse return")

    # Simple print of information
    def print(self):
        print(f"Moleculetype: {self.moleculetype_line}")

    def _pass(self, line):
        pass

    # For now, just store the line information
    def parser_moleculetype(self, line):
        print(f"molecule line {line}")
        self.moleculetype_line = line.rstrip()


# Create a chargeless topology
def create_chargeless_topology(path):
    r""" Create a charge-free topology (ITP files) for use in breaking up electrostatic forces
    """

    # Create the new directory
    Path(os.path.join(path, "toppar_noq")).mkdir(parents=True, exist_ok=True)

    # Read the various files from the original location with our ITP parser, then we
    # can make changes as we see fit
    for filename in os.listdir(os.path.join(path, "toppar")):
        #itp_parser = dfruit_ITPParser(os.path.join(path, "toppar", filename), True)
        #itp_parser.parse()
        #itp_parser.print()
        filestem = filename.split('.')[0]
        orifile = os.path.join(path, "toppar", filestem + ".itp")
        newfile = os.path.join(path, "toppar_noq", filestem + "_noq.itp")
        with open(orifile, 'r') as file1, open(newfile, 'w') as file2:
            section = "none"
            for line in file1:
                line_to_print = line.rstrip() # Line that we will print at the end

                if '[' in line and ']' in line:
                    section = line.split('[')[1].split(']')[0].strip()

                # We are in a known section, so process it
                else:
                    if section == "atoms":
                        # Here is where we change charges
                        if line.startswith(";") or line.startswith("\n"):
                            line_to_print = line.rstrip()
                        else:
                            tokens = line.rstrip().split()
                            # Change the charge
                            tokens[6] = 0.0
                            line_to_print = ""
                            for token in tokens:
                                line_to_print += "\t" + str(token)

                    elif section == "atomtypes":
                        # Charge is in a different location
                        if line.startswith(";") or line.startswith("\n"):
                            line_to_print = line.rstrip()
                        else:
                            tokens = line.rstrip().split()
                            # Change the charge
                            tokens[3] = 0.0
                            line_to_print = ""
                            for token in tokens:
                                line_to_print += "\t" + str(token)

                line_to_print += "\n"
                file2.write(line_to_print)

    # Finally, grab the topology file and write a new one (topol_noq.top)
    topol_ori_file = os.path.join(path, "topol.top")
    topol_noq_file = os.path.join(path, "topol_noq.top")
    with open(topol_ori_file, "r") as file1, open(topol_noq_file, "w") as file2:
        for line in file1:
            line_to_print = ""
            if line.startswith("#"):
                m0 = re.sub(r'(.itp)', r'_noq.itp', line.rstrip())
                m1 = re.sub(r'(toppar)', r'toppar_noq', m0)
                line_to_print = m1.rstrip() + "\n"
            else:
                line_to_print = line.rstrip() + "\n"
            
            file2.write(line_to_print)


