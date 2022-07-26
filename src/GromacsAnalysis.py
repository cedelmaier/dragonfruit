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

# Membrane curvature module is special
from membrane_curvature.base import MembraneCurvature

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from stylelib.common_styles import septin_runs_stl

def plot_contour_special(results, com, unit_cell, nbins, levels_, lmin, lmax, label, cm, savename):
    """
    """
    plt.style.use(septin_runs_stl)
    zoom_factor = 3

    fig, ax = plt.subplots(1, 1, figsize=(4, 3.5), dpi=200)
    max_ = np.max(results)
    print("Max_: {}".format(max_))
    rs = ndimage.zoom(results, zoom_factor, mode='wrap', order=1)
    levs = np.linspace(int(lmin), round(lmax), levels_)
    im = ax.contourf(rs, cmap=cm, origin='lower', levels=levs, alpha = 0.95, vmin=int(lmin), vmax=round(lmax))
    tcs = [lmin, lmax]

    scaled_coord = com / unit_cell[0:3] * nbins * zoom_factor
    ax.scatter(scaled_coord[1], scaled_coord[0], s=35, c='k', alpha=0.5)

    ax.set_aspect('equal')
    ax.axis('off')
    cbar = plt.colorbar(im, ticks=tcs, orientation='horizontal', ax=ax, shrink=0.7, aspect=10, pad=0.05)
    cbar.ax.tick_params(labelsize=4, width=0.5)
    cbar.set_label(label, fontsize=6, labelpad=2)

    fig.tight_layout()
    fig.savefig(savename, dpi=fig.dpi)

    return

def plot_contours_com(results, com, unit_cell, nbins, label, levels_, cm, savename):
    """
    Function used to plot contours of MembraneCurvature results.
    User can determine number of contour lines / regions (levels),
    label of the colorbar (label) and colormap (cmap).

    Parameters
    ----------
    results: list
        List with results by leaflets as elements [lower_leaflet, upper_leaflet]
    com:
        Center of mass of the protein
    unit_cell:
        Unit cell variable to plot
    label: str
        Label to add to colorbar.
    levels: int
        Determines number of contour lines.
    cmap: str
        Colormap to use in plot.

    """
    leaflets = ['Lower', 'Upper']

    zoom_factor = 3

    fig, [ax1, ax2] = plt.subplots(ncols=2, figsize=(4,3.5), dpi=200)
    max_ = np.max(results)
    for ax, rs, lf in zip((ax1, ax2), results, leaflets):
        rs = ndimage.zoom(rs, zoom_factor, mode='wrap', order=1)
        if np.min(rs) < 0 < np.max(rs):
            levs = np.linspace(-max_, max_, levels_)
            im = ax.contourf(rs, cmap=cm, origin='lower', levels=levs, alpha=0.95, vmin=-max_, vmax=max_)
            tcs = [-max_, 0, max_]
        else:
            levs = np.linspace(int(np.min(rs)), round(np.max(rs)), levels_)
            im = ax.contourf(rs, cmap=cm, origin='lower', levels=levs, alpha=0.95,  vmin=int(np.min(rs)), vmax=round(np.max(rs)))
            tcs = [int(np.min(rs)), round(np.max(rs))]

        # Get the location of the protein too
        # Also, need to reverse X-Y because of image conventions
        print("Xlim:               : {}".format(ax1.get_xlim()))
        print("Ylim:               : {}".format(ax1.get_ylim()))
        #scaled_coord = com / nbins * zoom_factor
        scaled_coord = com / unit_cell[0:3] * nbins * zoom_factor
        print("Center of mass (raw): {}".format(com))
        print("Unit cell           : {}".format(unit_cell[0:3]))
        print("Nbins:              : {}".format(nbins))
        print("Scaled coord (zoom) : {}".format(scaled_coord))
        ax.scatter(scaled_coord[1], scaled_coord[0], s=5, c='k', alpha=0.5)

        ax.set_aspect('equal')
        ax.set_title('{} Leaflet'.format(lf), fontsize=6)
        ax.axis('off')
        cbar = plt.colorbar(im, ticks=tcs, orientation='horizontal', ax=ax, shrink=0.7, aspect=10, pad=0.05)
        cbar.ax.tick_params(labelsize=4, width=0.5)
        cbar.set_label(label, fontsize=6, labelpad=2)

    # Save this off before returning
    fig.tight_layout()
    fig.savefig(savename, dpi=fig.dpi)

    return

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
        #self.filename_gromacs = os.path.basename(self.filename_structure).split('.')[0] + '.gro'

        #self.AnalyzeTrajectory()
        self.AnalyzeCurvature()
        #self.AnalyzeHelix()
        #self.AnalyzeDSSP()

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

    def AnalyzeCurvature(self):
        r""" Analyze the curvature of the membrane
        """
        print("GromacsAnalysis::AnalyzeCurvature")
        universe = mda.Universe(self.filename_structure, self.filename_trajectory)

        protein = universe.select_atoms("protein")
        print("\nProtein has {} residues." \
            "".format(protein.n_residues))

        lipids  = universe.select_atoms("not protein")
        print("\nLipids molecules include {} residues and {} atoms." \
            "".format(lipids.n_residues, lipids.n_atoms))

        P_headgroups = universe.select_atoms('name P')

        L_ah = leaf.LeafletFinder(universe, 'name P', cutoff = 20)
        ah_upper_leaflet = L_ah.groups(0)
        ah_lower_leaflet = L_ah.groups(1)

        leaflets = ['Lower', 'Upper']

        ah_upper_leaflet_P = ah_upper_leaflet.select_atoms('name P')
        ah_lower_leaflet_P = ah_lower_leaflet.select_atoms('name P')

        for name, new_lf in zip(leaflets, [ah_lower_leaflet_P, ah_upper_leaflet_P]):
            print("{} leaflet includes {} elements.".format(name, len(new_lf)))

        sel_upper = " ".join([str(r) for r in ah_upper_leaflet.residues.resids])
        sel_lower = " ".join([str(r) for r in ah_lower_leaflet.residues.resids])

        upper_string = "resid {} and name P".format(sel_upper)
        lower_string = "resid {} and name P".format(sel_lower)

        # Determine the unit_cell dimensions before doing the curvature calculation
        unit_cell = []
        for ts in universe.trajectory:
            unit_cell.append(universe.dimensions)
        unit_cell_new = np.stack(unit_cell)
        mean_unit_cell = np.mean(unit_cell_new, axis=0)

        # Set the x and y dimension sizes to the average for computation
        x_dim = (0, mean_unit_cell[0])
        y_dim = (0, mean_unit_cell[1])

        # Do the curvature analysis
        nbins = 6
        curvature_upper_leaflet = MembraneCurvature(universe,
                                                    select = upper_string,
                                                    n_x_bins = nbins,
                                                    n_y_bins = nbins,
                                                    x_range = x_dim,
                                                    y_range = y_dim,
                                                    wrap = True).run()
        curvature_lower_leaflet = MembraneCurvature(universe,
                                                    select = lower_string,
                                                    n_x_bins = nbins,
                                                    n_y_bins = nbins,
                                                    x_range = x_dim,
                                                    y_range = y_dim,
                                                    wrap = True).run()
        helix_com = []
        for ts in universe.trajectory:
            helix_com.append(protein.center_of_mass())

        helix_com_new = np.stack(helix_com)

        mean_helix_com = np.mean(helix_com_new, axis=0)

        # Do the plotting here until we can clean up this mess
        avg_surface_upper_leaflet = curvature_upper_leaflet.results.average_z_surface
        avg_surface_lower_leaflet = curvature_lower_leaflet.results.average_z_surface

        avg_mean_upper_leaflet = curvature_upper_leaflet.results.average_mean
        avg_mean_lower_leaflet = curvature_lower_leaflet.results.average_mean

        avg_gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian
        avg_gaussian_lower_leaflet = curvature_lower_leaflet.results.average_gaussian

        # Process everything different
        plot_contour_special(avg_mean_lower_leaflet, mean_helix_com, mean_unit_cell, nbins, 30, -1.0, 1.0, "$H$ (Å$^{-1}$)", "bwr", 'avg_mean_lower.pdf')
        plot_contour_special(avg_mean_upper_leaflet, mean_helix_com, mean_unit_cell, nbins, 30, -1.0, 1.0, "$H$ (Å$^{-1}$)", "bwr", 'avg_mean_upper.pdf')
        plot_contour_special(avg_gaussian_lower_leaflet, mean_helix_com, mean_unit_cell, nbins, 35, -2.0, 2.0, "$K$ (Å$^{-2}$)", "PiYG", 'avg_gaussian_lower.pdf')
        plot_contour_special(avg_gaussian_upper_leaflet, mean_helix_com, mean_unit_cell, nbins, 35, -2.0, 2.0, "$K$ (Å$^{-2}$)", "PiYG", 'avg_gaussian_upper.pdf')

        print("Exiting here!")
        sys.exit(1)

        ## Try to calculate the principle curvatures from this
        #from numpy.polynomial import Polynomial
        #avg_K1_upper = avg_surface_upper_leaflet
        #for idx, _ in np.ndenumerate(avg_surface_upper_leaflet):
        #    array_vals = np.array([avg_gaussian_upper_leaflet[idx], -2.0*avg_mean_upper_leaflet[idx], 1.0])
        #    print("array vals:      {}".format(array_vals))
        #    pval = Polynomial(array_vals)
        #    print("roots:           {}".format(pval.roots()))

        #    sys.exit(1)
        

        plot_contours_com([avg_surface_lower_leaflet, avg_surface_upper_leaflet], mean_helix_com, mean_unit_cell, nbins, 'Surface', 35, 'YlGnBu_r', 'surface_height.pdf')
        #plt.show()

        plot_contours_com([avg_mean_lower_leaflet, avg_mean_upper_leaflet], mean_helix_com, mean_unit_cell, nbins, "$H$ (Å$^{-1}$)", 30, "bwr", 'mean_curvature.pdf')
        #plt.show()

        plot_contours_com([avg_gaussian_lower_leaflet, avg_gaussian_upper_leaflet], mean_helix_com, mean_unit_cell, nbins, "$K$ (Å$^{-2}$)", 35, "PiYG", 'gaussian_curvature.pdf')
        #plt.show()

        print("Premature exit")
        sys.exit(1)


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

        ## Plot the helicity of the protein
        #fig, axarr = plt.subplots(1, 1, figsize = (15, 20))
        #self.graph_helicity(axarr)
        #fig.tight_layout()
        #fig.savefig('gromacs_helicity.pdf', dpi=fig.dpi)

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
