#!/usr/bin/env python3

# XXX: Put a license here

"""Main analysis script for membranes with AH domains"""

import argparse
import pickle
import os
import sys

import MDAnalysis as mda
import MDAnalysis.transformations as trans
from MDAnalysis.analysis import helix_analysis as hel

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from stylelib.common_styles import septin_runs_stl

def parse_args():
    parser = argparse.ArgumentParser(prog='GromacsAnalysis.py')

    # General options
    parser.add_argument('-tpr', '--structure', required=True, type=str,
            help = 'Gromacs TPR topology/structure file')
    parser.add_argument('-xtc', '--trajectory', required=True, type=str,
            help = 'Gromacs TPR trajectory file')
    parser.add_argument('-d', '--workdir', action = 'store_true',
            help = 'Working directory')


    opts = parser.parse_args()
    return opts

class GromacsAnalysis(object):
    r""" Gromacs analysis
    """
    def __init__(self, opts):
        self.opts = opts
        self.cwd = os.getcwd()

        self.ReadOpts()

        self.ProgOpts()

    def ReadOpts(self):
        if not self.opts.workdir:
            self.opts.workdir = os.path.abspath(self.cwd)
        elif not os.path.exists(self.opts.workdir):
            raise IOError("Working directory {} does not exist.".format(self.opts.workdir))
        else:
            self.opts.workdir = os.path.abspath(self.opts.workdir)

    def ProgOpts(self):
        r""" Run selected commands
        """

        self.Analyze()

    def Analyze(self):
        r""" Run analysis on a single simulation
        """
        self.path = os.path.abspath(self.opts.workdir)
        self.name = self.path.split('/')[-1]

        self.filename_structure = os.path.join(self.path, self.opts.structure)
        self.filename_trajectory = os.path.join(self.path, self.opts.trajectory)

        self.AnalyzeTrajectory()
        self.AnalyzeHelix()

        self.Graph()
        self.WriteData()

    def AnalyzeTrajectory(self):
        r""" Analyze trajectory information, unwrapped coordinates
        """
        print("GromacsAnalysis::AnalyzeTrajectory")
        # Create the universe
        traj_universe = mda.Universe(self.filename_structure, self.filename_trajectory)

        dopc_atoms = traj_universe.select_atoms('resname DOPC')
        plpi_atoms = traj_universe.select_atoms('resname PLPI')
        lipid_atoms = traj_universe.select_atoms('resname DOPC or resname PLPI')
        helix_atoms = traj_universe.select_atoms('protein')
        not_protein = traj_universe.select_atoms('not protein')
        not_solvent = traj_universe.select_atoms('resname DOPC or resname PLPI or protein')
        solvent = traj_universe.select_atoms('not (resname DOPC or resname PLPI or protein)')

        # Wrap everything into a periodic box
        #transforms = [trans.unwrap(not_solvent),
        #              trans.center_in_box(not_solvent, wrap=True),
        #              trans.wrap(solvent)]
        transforms = [trans.unwrap(not_solvent)]
        traj_universe.trajectory.add_transformations(*transforms)
        
        # Do the analysis on the trajectory
        times = []
        dopc_com = []
        plpi_com = []
        lipid_com = []
        helix_com = []
        for ts in traj_universe.trajectory:
            times.append(traj_universe.trajectory.time)
            dopc_com.append(dopc_atoms.center_of_mass())
            plpi_com.append(plpi_atoms.center_of_mass())
            lipid_com.append(lipid_atoms.center_of_mass())
            helix_com.append(helix_atoms.center_of_mass())

        # Split up into a pandas dataframe for viewing
        dopc_com_df = pd.DataFrame(dopc_com, columns = ['dopc_x', 'dopc_y', 'dopc_z'], index=times)
        dopc_com_df.index.name = 'Time(ps)'
        plpi_com_df = pd.DataFrame(plpi_com, columns = ['plpi_x', 'plpi_y', 'plpi_z'], index=times)
        plpi_com_df.index.name = 'Time(ps)'
        lipid_com_df = pd.DataFrame(lipid_com, columns = ['lipid_x', 'lipid_y', 'lipid_z'], index=times)
        lipid_com_df.index.name = 'Time(ps)'
        helix_com_df = pd.DataFrame(helix_com, columns = ['helix_x', 'helix_y', 'helix_z'], index=times)
        helix_com_df.index.name = 'Time(ps)'
        
        dfs = []
        dfs.append(dopc_com_df)
        dfs.append(plpi_com_df)
        dfs.append(lipid_com_df)
        dfs.append(helix_com_df)
        
        self.master_df = pd.concat(dfs, axis=1)

    def AnalyzeHelix(self):
        r""" Analyze helix parameters
        """
        print(f"GromacsAnalysis::AnalyzeHelix")
        # Create a new universe for coordinate transfomrs
        traj_universe_helix = mda.Universe(self.filename_structure, self.filename_trajectory)

        # Do the transformation for the COM coordinates of the helix to stay in the box
        all_protein = traj_universe_helix.select_atoms('protein')
        not_protein = traj_universe_helix.select_atoms('not protein')

        # Wrap everything into a periodic box
        transforms = [trans.unwrap(all_protein),
                      trans.center_in_box(all_protein, wrap=True),
                      trans.wrap(not_protein)]

        traj_universe_helix.trajectory.add_transformations(*transforms)

        # Run the helix analysis
        self.helix_analysis = hel.HELANAL(traj_universe_helix, select='name CA and resid 1-18').run()

    def Graph(self):
        r""" Graph results
        """
        print(f"GromacsAnalysis::Graph")
        plt.style.use(septin_runs_stl)

        # Plot the Z distance versus time
        fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
        self.graph_zdist_com(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_zdist.pdf', dpi=fig.dpi)

        # Plot helix information
        fig, axarr = plt.subplots(2, 1, figsize = (15, 10))
        self.graph_helix_analysis(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_helix_analysis.pdf', dpi=fig.dpi)

    def WriteData(self):
        r""" Write the data into either CSV files for pandas or pickles for internal objects
        """
        # Dump the CSV files
        hd5_filename = os.path.join(self.path, self.name + '.h5')
        self.master_df.to_hdf(hd5_filename, key='master_df', mode='w')

        # Dump the pickle files
        pkl_filename = os.path.join(self.path, self.name + '.pickle')
        with open(pkl_filename, 'wb') as f:
            pickle.dump(self.helix_analysis,f)


    def graph_zdist_com(self, axarr):
        r""" Plot the Z distance between the lipid and helix COM
        """
        z = np.abs(self.master_df['helix_z'] - self.master_df['lipid_z'])
        #z.plot()
        axarr.plot(z)
        axarr.set_xlabel('Frame')
        axarr.set_ylabel('Z distance (Angstroms)')

    def graph_helix_analysis(self, axarr):
        r""" Plot results of the helix analysis
        """
        axarr[0].plot(self.helix_analysis.results.local_twists.mean(axis=1))
        axarr[0].set_xlabel('Frame')
        axarr[0].set_ylabel('Average twist (degrees)')

        axarr[1].plot(self.helix_analysis.results.local_nres_per_turn.mean(axis=1))
        axarr[1].set_xlabel('Frame')
        axarr[1].set_ylabel('Average residues per turn')





##########################################
if __name__ == "__main__":
    opts = parse_args()
    x = GromacsAnalysis(opts)
