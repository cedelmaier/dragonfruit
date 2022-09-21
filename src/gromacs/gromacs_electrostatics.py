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
from dragonfruit_io import create_chargeless_topology, read_xvg

class GromacsElectrostaticsAnalyzer(object):
    def __init__(self, path, opts):
        self.verbose = opts.verbose
        self.trace = opts.trace
        if self.verbose: print("GromacsElectrostaticsAnalyzer::__init__")

        self.path = os.path.abspath(path)

        self.electrostatic_dfs = {}

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

    def GenerateForces(self, force_type, trajectory_file):
        r""" Generate the forces according to the decomposition MDP file
        """
        if self.verbose: print(f"GromacsElectrostaticsAnalyzer::GenerateForces {force_type}")

        # Remember, everything is under the auspices of decomp_base.mdp
        filename_base = "decomp_" + force_type

        # Remove previous files that were generated
        tpr_file = filename_base + ".tpr"
        trr_file = filename_base + ".trr"
        edr_file = filename_base + ".edr"
        xvg_file = filename_base + ".xvg"
        log_file = filename_base + ".log"

        protein_force_file = "protein_force_" + force_type + ".xvg"
        interaction_energy_file = "interaction_energy_" + force_type + ".xvg"
        new_traj_filename = trajectory_file.split(".")[0] + "_resampled10.xtc"

        if os.path.exists(os.path.join(self.path, tpr_file)): os.remove(os.path.join(self.path, tpr_file))
        if os.path.exists(os.path.join(self.path, trr_file)): os.remove(os.path.join(self.path, trr_file))
        if os.path.exists(os.path.join(self.path, edr_file)): os.remove(os.path.join(self.path, edr_file))
        if os.path.exists(os.path.join(self.path, xvg_file)): os.remove(os.path.join(self.path, xvg_file))
        if os.path.exists(os.path.join(self.path, log_file)): os.remove(os.path.join(self.path, log_file))

        if os.path.exists(os.path.join(self.path, protein_force_file)): os.remove(os.path.join(self.path, protein_force_file))
        if os.path.exists(os.path.join(self.path, interaction_energy_file)): os.remove(os.path.join(self.path, interaction_energy_file))
        if os.path.exists(os.path.join(self.path, new_traj_filename)): os.remove(os.path.join(self.path, new_traj_filename))

        topo_file = ""
        if force_type == "all":
            topo_file = "topol.top"
        elif force_type == "noq":
            topo_file = "topol_noq.top"
        else:
            print(f"ERROR: force_type {force_type} requested, not implemented, exiting!")
            sys.exit(1)

        # Reduce the size of the trajectory file, as it gets out of hand quickly
        gmx_trjconv_p = subprocess.Popen(["gmx", "trjconv", "-s", "step6.6_equilibration.tpr", "-f", trajectory_file, "-o", new_traj_filename, "-n", "density_groups.ndx", "-skip", "10"], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        gmx_trjconv_p.stdin.write(b'0\n')
        trjconv_output, trjconv_errors = gmx_trjconv_p.communicate()
        gmx_trjconv_p.wait()
        if self.trace: print(trjconv_errors.decode("utf-8"))
        if self.trace: print(trjconv_output.decode("utf-8"))

        # Create the run file
        gmx_grompp_p = subprocess.Popen(["gmx", "grompp", "-f", "decomp_base.mdp", "-c", "step6.6_equilibration.gro", "-t", "step6.6_equilibration.cpt", "-p", topo_file, "-n", "density_groups.ndx", "-o", tpr_file, "-maxwarn", "1"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        grompp_output, grompp_errors = gmx_grompp_p.communicate()
        gmx_grompp_p.wait()
        if self.trace: print(grompp_errors.decode("utf-8"))
        if self.trace: print(grompp_output.decode("utf-8"))

        # Execute force extract in single threaded mode
        gmx_mdrun_p = subprocess.Popen(["gmx", "mdrun", "-deffnm", filename_base, "-rerun", new_traj_filename], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        mdrun_output, mdrun_errors = gmx_mdrun_p.communicate()
        gmx_mdrun_p.wait()
        if self.trace: print(mdrun_errors.decode("utf-8"))
        if self.trace: print(mdrun_output.decode("utf-8"))

        # Finally, can generate the forces
        gmx_traj_p = subprocess.Popen(["gmx", "traj", "-f", trr_file, "-s", tpr_file, "-n", "density_groups.ndx", "-of", protein_force_file], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        gmx_traj_p.stdin.write(b'1\n')
        traj_output, traj_errors = gmx_traj_p.communicate()
        gmx_traj_p.wait()
        if self.trace: print(traj_errors.decode("utf-8"))
        if self.trace: print(traj_output.decode("utf-8"))

        # Do the energies too, while we're at it
        gmx_energy_p = subprocess.Popen(["gmx", "energy", "-f", filename_base, "-o", interaction_energy_file], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        gmx_energy_p.stdin.write(b'21\n22\n23\n24\n\n')
        energy_output, energy_errors = gmx_energy_p.communicate()
        gmx_energy_p.wait()
        if self.trace: print(energy_errors.decode("utf-8"))
        if self.trace: print(energy_output.decode("utf-8"))

        # Now, read in the XVG file for the force, and load it into a dataframe
        metadata, num_data, df = read_xvg(protein_force_file, force_type)
        self.electrostatic_dfs[force_type] = df

        # Remove all of the files at the end of the run to save space, we have everything we want in the dataframe for now
        if os.path.exists(os.path.join(self.path, tpr_file)): os.remove(os.path.join(self.path, tpr_file))
        if os.path.exists(os.path.join(self.path, trr_file)): os.remove(os.path.join(self.path, trr_file))
        if os.path.exists(os.path.join(self.path, edr_file)): os.remove(os.path.join(self.path, edr_file))
        if os.path.exists(os.path.join(self.path, xvg_file)): os.remove(os.path.join(self.path, xvg_file))
        if os.path.exists(os.path.join(self.path, log_file)): os.remove(os.path.join(self.path, log_file))

        if os.path.exists(os.path.join(self.path, protein_force_file)): os.remove(os.path.join(self.path, protein_force_file))
        if os.path.exists(os.path.join(self.path, interaction_energy_file)): os.remove(os.path.join(self.path, interaction_energy_file))
        if os.path.exists(os.path.join(self.path, new_traj_filename)): os.remove(os.path.join(self.path, new_traj_filename))

        if self.verbose: print(f"GromacsElectrostaticsAnalyzer::GenerateForces {force_type} return")

    def CalculateElectrostaticForces(self):
        r""" Calculate the electrostatic forces given the all and noq force dataframes
        
        Requires that both the all and noq force dataframes exist
        """
        if self.verbose: print("GromacsElectrostaticsAnalyzer::CalculateElectrostaticForces")
        # Create the q dataframe
        df_q = self.electrostatic_dfs["all"] - self.electrostatic_dfs["noq"].values
        # Rename the columns in this last one
        df_q.rename(columns = lambda x: x.replace("all", "q"), inplace = True)
        self.electrostatic_df = pd.concat([self.electrostatic_dfs["all"], self.electrostatic_dfs["noq"], df_q], axis = 1)

        if self.verbose: print("GromacsElectrostaticsAnalyzer::CalculateElectrostaticForces return")

    def WriteDecompositionMDP(self):
        r""" Write decomposition MDP file for gromacs
        """
        if self.verbose: print("GromacsElectrostaticsAnalyzer::WriteDecompositionMDP")
        # Create the decomposition MDP file
        with open(os.path.join(self.path, "decomp_base.mdp"), "w") as fhandle:
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
        if os.path.exists(os.path.join(self.path, "density_groups.ndx")): os.remove(os.path.join(self.path, "density_groups.ndx"))
        with open(os.path.join(self.path, "create_density_map.txt"), "r") as in_stream:
                gmx_density_p = subprocess.Popen(["gmx", "make_ndx", "-f", "step7_20.tpr", "-o", "density_groups.ndx"], stdin = in_stream, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                g_output, g_errors = gmx_density_p.communicate()
                gmx_density_p.wait()
                if self.trace: print(g_errors.decode("utf-8"))
                if self.trace: print(g_output.decode("utf-8"))

        if self.verbose: print("GromacsElectrostaticsAnalyzer::WriteCreateDensityMapTXT return")

