#!/usr/bin/env python3

# XXX: Put a license here

"""Main analysis script for membranes with AH domains"""

import argparse
import pickle
import os
import sys
import time

# MD Analysis is used for leaflet identification, transformations
import MDAnalysis as mda
import MDAnalysis.transformations as trans
from MDAnalysis.analysis import helix_analysis as hel
from MDAnalysis.analysis import leaflet as leaf

# MDTraj has a native DSSP module
import mdtraj as md

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from stylelib.common_styles import septin_runs_stl

def parse_args():
    parser = argparse.ArgumentParser(prog='GromacsAnalysis.py')

    # General options
    parser.add_argument('-top', '--structure', required=True, type=str,
            help = 'Gromacs topology/structure file')
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
        self.start_time = time.time()
        self.opts = opts
        self.cwd = os.getcwd()

        self.ReadOpts()

        self.ProgOpts()

        print("--- %s seconds ---" % (time.time() - self.start_time))

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

        # Create a reader depending on what we are reading in
        self.filename_structure = os.path.join(self.path, self.opts.structure)
        self.filename_trajectory = os.path.join(self.path, self.opts.trajectory)
        self.filename_gromacs = os.path.basename(self.filename_structure).split('.')[0] + '.gro'

        self.AnalyzeTrajectory()
        #self.AnalyzeHelix()
        self.AnalyzeDSSP()

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

        # Try to get the head groups of the two leaflets for analysis too
        # Has an implicit cutoff at 15.0
        L = leaf.LeafletFinder(traj_universe, 'name P*')
        leaflet0 = L.groups(0)
        leaflet1 = L.groups(1)

        # Unwrap/wrap the protein so that it doesn't have issues with the PBC
        #transforms = [trans.unwrap(helix_atoms)]
        transforms = [trans.unwrap(helix_atoms),
                      trans.unwrap(lipid_atoms),
                      trans.center_in_box(lipid_atoms, wrap=True),
                      trans.wrap(solvent),
                      trans.wrap(lipid_atoms)]
        traj_universe.trajectory.add_transformations(*transforms)

        # Wrap everything into a periodic box
        #transforms = [trans.unwrap(not_solvent),
        #              trans.center_in_box(not_solvent, wrap=True),
        #              trans.wrap(solvent)]
        #transforms = [trans.unwrap(not_solvent)]
        #traj_universe.trajectory.add_transformations(*transforms)
        
        # Do the analysis on the trajectory
        times = []
        dopc_com = []
        plpi_com = []
        lipid_com = []
        helix_com = []
        leaflet0_com = []
        leaflet1_com = []
        unit_cell = []
        for ts in traj_universe.trajectory:
            times.append(traj_universe.trajectory.time)
            dopc_com.append(dopc_atoms.center_of_mass())
            plpi_com.append(plpi_atoms.center_of_mass())
            lipid_com.append(lipid_atoms.center_of_mass())
            helix_com.append(helix_atoms.center_of_mass())
            leaflet0_com.append(leaflet0.center_of_mass())
            leaflet1_com.append(leaflet1.center_of_mass())

            # Get the unit cell as well
            unit_cell.append(traj_universe.dimensions)
            print(traj_universe.dimensions)
            sys.exit(1)

        # Save off the times for other uses!
        self.times = times

        # Split up into a pandas dataframe for viewing
        dopc_com_df = pd.DataFrame(dopc_com, columns = ['dopc_x', 'dopc_y', 'dopc_z'], index=times)
        dopc_com_df.index.name = 'Time(ps)'
        plpi_com_df = pd.DataFrame(plpi_com, columns = ['plpi_x', 'plpi_y', 'plpi_z'], index=times)
        plpi_com_df.index.name = 'Time(ps)'
        lipid_com_df = pd.DataFrame(lipid_com, columns = ['lipid_x', 'lipid_y', 'lipid_z'], index=times)
        lipid_com_df.index.name = 'Time(ps)'
        helix_com_df = pd.DataFrame(helix_com, columns = ['helix_x', 'helix_y', 'helix_z'], index=times)
        helix_com_df.index.name = 'Time(ps)'
        leaflet0_com_df = pd.DataFrame(leaflet0_com, columns = ['leaflet0_x', 'leaflet0_y', 'leaflet0_z'], index=times)
        leaflet0_com_df.index.name = 'Time(ps)'
        leaflet1_com_df = pd.DataFrame(leaflet1_com, columns = ['leaflet1_x', 'leaflet1_y', 'leaflet1_z'], index=times)
        leaflet1_com_df.index.name = 'Time(ps)'
        unit_cell_df = pd.DataFrame(unit_cell, columns = ['unit_cell_x', 'unit_cell_y', 'unit_cell_z', 'unit_cell_alpha', 'unit_cell_beta', 'unit_cell_gamma'], index = times)
        unit_cell_df.index.name = 'Time(ps)'

        dfs = []
        dfs.append(dopc_com_df)
        dfs.append(plpi_com_df)
        dfs.append(lipid_com_df)
        dfs.append(helix_com_df)
        dfs.append(leaflet0_com_df)
        dfs.append(leaflet1_com_df)
        dfs.append(unit_cell_df)
        
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

    def AnalyzeDSSP(self):
        r""" Analyze DSSP parameters for secondary structure
        """
        print(f"GromacsAnalysis::AnalyzeDSSP")
        # Iteratively load as otherwise this takes forever!
        current_chunk = 0
        # Create an array to store the helicity results
        helicity_results = []
        for chunk in md.iterload(self.filename_trajectory, top=self.filename_gromacs):
            ## Easy exit to not process the whole thing
            #if current_chunk >=1:
            #    break

            topology = chunk.topology
            helix_selection = topology.select('protein')
            helix_traj = chunk.atom_slice(helix_selection)

            dssp_results = md.compute_dssp(helix_traj)

            # Define a helicity for every frame based on this
            for iframe in dssp_results:
                helicity_results.append(iframe)

            current_chunk += 1

        helicity_results = np.array(helicity_results)

        self.helicity = np.zeros(helicity_results.shape[0])
        for iframe,val in enumerate(helicity_results):
            # Get the helicity of the combined states
            h_count = (val == 'H').sum() + (val == 'G').sum() + (val == 'I').sum()
            h_total = val.shape

            current_helicity = h_count / h_total
            self.helicity[iframe] = current_helicity[0]

        # Add the helicity to the master DF
        helicity_df = pd.DataFrame(self.helicity, columns = ['helicity'], index=self.times)
        helicity_df.index.name = 'Time(ps)'

        self.master_df = pd.concat([self.master_df, helicity_df], axis=1)

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

        # Try doing a plot with the lipid head groups included
        fig, axarr = plt.subplots(1, 1, figsize=(15,10))
        self.graph_zpos_com(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_zpos.pdf', dpi=fig.dpi)

        # Plot the helicity of the protein
        fig, axarr = plt.subplots(1, 1, figsize = (15, 20))
        self.graph_helicity(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_helicity.pdf', dpi=fig.dpi)

        ## Plot helix information
        #fig, axarr = plt.subplots(2, 1, figsize = (15, 10))
        #self.graph_helix_analysis(axarr)
        #fig.tight_layout()
        #fig.savefig('gromacs_helix_analysis.pdf', dpi=fig.dpi)

    def WriteData(self):
        r""" Write the data into either CSV files for pandas or pickles for internal objects
        """
        # Dump the CSV files
        hd5_filename = os.path.join(self.path, self.name + '.h5')
        self.master_df.to_hdf(hd5_filename, key='master_df', mode='w')

        ## Dump the pickle files
        #pkl_filename = os.path.join(self.path, self.name + '.pickle')
        #with open(pkl_filename, 'wb') as f:
        #    pickle.dump(self.helix_analysis,f)


    def graph_zdist_com(self, axarr):
        r""" Plot the Z distance between the lipid and helix COM
        """
        z = np.abs(self.master_df['helix_z'] - self.master_df['lipid_z'])
        axarr.plot(z)
        axarr.set_xlabel('Frame')
        axarr.set_ylabel('Z distance (Angstroms)')

    def graph_zpos_com(self, axarr):
        r""" Plot the Z positino of the lipid head groups and the protein
        """
        z_protein = self.master_df['helix_z']
        z_leaf0 = self.master_df['leaflet0_z']
        z_leaf1 = self.master_df['leaflet1_z']
        z_lipid = self.master_df['lipid_z']
        # Subtract off the position of the lipid COM from everybody else
        z_protein = z_protein - z_lipid
        z_leaf0 = z_leaf0 - z_lipid
        z_leaf1 = z_leaf1 - z_lipid
        axarr.plot(z_protein, color = 'b')
        axarr.plot(z_leaf0, color = 'k')
        axarr.plot(z_leaf1, color = 'k')
        axarr.set_xlabel('Frame')
        axarr.set_ylabel('Z position (Angstroms)')

    def graph_helicity(self, axarr):
        r""" Plot the helicity (0 to 1) as a function of time
        """
        helicity = self.master_df['helicity']
        axarr.plot(helicity, color = 'k')
        axarr.set_xlabel('Frame')
        axarr.set_ylabel('Helicity (AU)')
        axarr.set_ylim([0.0, 1.05])

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
