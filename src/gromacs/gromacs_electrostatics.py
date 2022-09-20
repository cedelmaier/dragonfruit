#!/usr/bin/env python3

# XXX: Put a license here

""" Electorstatics specific class for gromacs AH domains """

import os
import subprocess
import sys

import numpy as np
import pandas as pd

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from stylelib.common_styles import septin_runs_stl
from dragonfruit_io import create_chargeless_topology

class GromacsElectrostaticsAnalyzer(object):
    def __init__(self, path, opts):
        self.verbose = opts.verbose
        if self.verbose: print("GromacsElectrostaticsAnalyzer::__init__")
        self.path = os.path.abspath(path)
        if self.verbose: print("GromacsElectrostaticsAnalyzer::__init__ return")

    def CreateChargelessTopology(self, lipids):
        r""" Setup the charge-less electrostatics analysis, along with any changes to files needed

        Needs to know about the lipids used in the simulation, as well as the path
        """
        if self.verbose: print("GromacsElectrostaticsAnalyzer::CreateChargelessTopology")

        # Get what kind of lipids we are using
        self.lipids = lipids
        if self.verbose: print(f"  GromacsElectrostaticsAnalyzer lipid types: {self.lipids}")

        # Create the chargless topology directory
        create_chargeless_topology(self.path)

        # Create the chargless decomposition files
        self.WriteDecompositionMDP()
        self.WriteCreateDensityMapTXT()

        if self.verbose: print("GromacsElectrostaticsAnalyzer::CreateChargelessTopology return")

    def GenerateForces(self, force_type):
        r""" Generate the forces according to the decomposition MDP file
        """
        if self.verbose: print("GromacsElectrostaticsAnalyzer::GenerateForces")
        if self.verbose: print("GromacsElectrostaticsAnalyzer::GenerateForces return")

    def WriteDecompositionMDP(self):
        r""" Write decomposition MDP file for gromacs
        """
        if self.verbose: print("GromacsElectrostaticsAnalyzer::WriteDecompositionMDP")
        # Create the decomposition MDP file
        with open(os.path.join(self.path, "decomp_all.mdp"), "w") as fhandle:
            fhandle.write(r'''integrator              = md
dt                      = 0.002
nsteps                  = 100000000
nstenergy               = 50000
nstlog                  = 50000
nstxout-compressed      = 50000
nstfout                 = 50000
energygrps              = Protein lipids
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = Nose-Hoover
tc_grps                 = System
tau_t                   = 1.0
ref_t                   = 300
;
pcoupl                  = Parrinello-Rahman
pcoupltype              = semiisotropic
tau_p                   = 5.0
compressibility         = 4.5e-5  4.5e-5
ref_p                   = 1.0     1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;''')
        if self.verbose: print("GromacsElectrostaticsAnalyzer::WriteDecompositionMDP return")


    def WriteCreateDensityMapTXT(self):
        r""" Write create_density_group.txt file for specific lipids
        """
        if self.verbose: print("GromacsElectrostaticsAnalyzer::WriteCreateDensityMapTXT")

        # Figure out what kind of lipids we have, and create the corresponding file
        if 'DOPC' in self.lipids:
            with open(os.path.join(self.path, "create_density_map.txt"), "w") as fhandle:
                fhandle.write(r'''r DOPC PLPI
name 18 lipids
1 | 18
name 19 helix_lipids
13 & a P | a O11 | a O12 | a O13 | a O14
name 20 dopc_phosphates
13 & a C21 | a O22 | a C31 | a O32
name 21 dopc_carbonyls
13 & a C218 | a C318
name 22 dopc_terminalmethyls
14 & a P | a O11 | a O12 | a O13 | a O14
name 23 plpi_phosphates
14 & a C21 | a O22 | a C31 | a O32
name 24 plpi_carbonyls
14 & a C218 | a C316
name 25 plpi_terminalmethyls
20 | 23
name 26 all_phosphates
21 | 24
name 27 all_carbonyls
22 | 25
name 28 all_terminalmethyls
q
''')
        else:
            print(f"ERROR: Currently do not have this type of lipid composition supported, exiting!")
            sys.exit(1)

        # Now actually create the density map in gromacs
        print(f"WARNING: Doing something very dirty here in electorstatics, check later!")
        with open(os.path.join(self.path, "create_density_map.txt"), "r") as in_stream:
                gmx_density_p = subprocess.Popen(["gmx", "make_ndx", "-f", "step7_20.tpr", "-o", "density_groups.ndx"], stdin = in_stream, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                g_output, g_errors = gmx_density_p.communicate()
                gmx_density_p.wait()
                if self.verbose: print(g_errors.decode("utf-8"))
                if self.verbose: print(g_output.decode("utf-8"))

        if self.verbose: print("GromacsElectrostaticsAnalyzer::WriteCreateDensityMapTXT return")

