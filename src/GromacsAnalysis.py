#!/usr/bin/env python3

# XXX: Put a license here

"""Main analysis script for membranes with AH domains"""

import argparse
import logging
import pickle
import os
import sys
import time

# MD Analysis is used for leaflet identification, transformations
import MDAnalysis as mda
import MDAnalysis.transformations as trans
from MDAnalysis.analysis import helix_analysis as hel
from MDAnalysis.analysis import leaflet as leaf
from MDAnalysis.lib.log import ProgressBar

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

        ## Get the location of the protein too
        ## Also, need to reverse X-Y because of image conventions
        #print("Xlim:               : {}".format(ax1.get_xlim()))
        #print("Ylim:               : {}".format(ax1.get_ylim()))
        ##scaled_coord = com / nbins * zoom_factor
        #scaled_coord = com / unit_cell[0:3] * nbins * zoom_factor
        #print("Center of mass (raw): {}".format(com))
        #print("Unit cell           : {}".format(unit_cell[0:3]))
        #print("Nbins:              : {}".format(nbins))
        #print("Scaled coord (zoom) : {}".format(scaled_coord))
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

def unpack_curvature(df, leaflet, measurement):
    r""" Unpack the curvature dataframe back into an array structure
    """
    # Unpack the curvature dataframe
    nframes = len(df.index)
    # Look for the average z surface (to get the bins)
    # XXX: This is kinda hacky for now, but it works....
    target_name = leaflet + '_' + measurement
    bin_result = list(filter(lambda x: x.startswith(target_name), df.columns))
    nbins_string = bin_result[-1]
    nbins = np.int32(nbins_string.split('_')[-1]) + 1

    # Get the filtered columns from the dataframe
    filtered_df = df.filter(regex = target_name + '*')

    # Get the reshaped version of the array
    curvature_measurement = filtered_df.to_numpy().reshape(nframes, nbins, nbins)

    return curvature_measurement

def parse_args():
    parser = argparse.ArgumentParser(prog='GromacsAnalysis.py')

    # General options
    parser.add_argument('-top', '--structure', required=True, type=str,
            help = 'Gromacs .tpr file')
    parser.add_argument('-xtc', '--trajectory', required=True, type=str,
            help = 'Gromacs .xtc file')
    parser.add_argument('-gro', '--gromacs', required=True, type=str,
            help = 'Gromacs .gro file')
    parser.add_argument('-d', '--workdir', action = 'store_true',
            help = 'Working directory')

    # Flow control options
    parser.add_argument('-A', '--analysis', action = 'store_true',
            help = 'Analyze data from MD simulation')
    parser.add_argument('-G', '--graph', action = 'store_true',
            help = 'Graph data from MD simulation')
    parser.add_argument('-W', '--write', action = 'store_true',
            help = 'Write data to HD5 and pickle files')
    parser.add_argument('-F', '--force', action = 'store_true',
            help = 'Force complete analysis of simulation(s)')

    # Add verbosity control
    parser.add_argument('-v', '--verbose', action="store_true",
            help = 'Verbose output')

    opts = parser.parse_args()
    return opts

class GromacsAnalysis(object):
    r""" Gromacs analysis
    """
    def __init__(self, opts):
        # Try to set up the log facilities
        mda_logger= logging.getLogger('MDAnalysis')
        if mda_logger.hasHandlers():
            mda_logger.handlers = []
        mda.start_logging()
        self.logger = logging.getLogger('MDAnalysis')

        self.opts = opts
        self.cwd = os.getcwd()

        self.ReadOpts()

        # Setup internal variables
        self.time_dfs = []
        self.lipid_set = {"DOPC",
                          "PLPI",
                          "CHL1",
                          "POPC",
                          "DOPE",
                          "DOPS",
                          "SAPI24",
                          "SAPI25"}

        self.ProgOpts()

    def ReadOpts(self):
        if not self.opts.workdir:
            self.opts.workdir = os.path.abspath(self.cwd)
        elif not os.path.exists(self.opts.workdir):
            raise IOError("Working directory {} does not exist.".format(self.opts.workdir))
        else:
            self.opts.workdir = os.path.abspath(self.opts.workdir)

        self.verbose = False
        if self.opts.verbose:
            self.verbose = True

        # Set up path information
        self.path = os.path.abspath(self.opts.workdir)
        self.name = self.path.split('/')[-1]
        if self.verbose: print(f"  Name: {self.name}")

        # Set up some common names
        self.hd5_name = self.name + '.h5'
        self.pickle_name = self.name + '.pickle'

    def ProgOpts(self):
        r""" Run selected commands
        """

        if self.opts.analysis:
            self.RunAnalysis()
        if self.opts.graph:
            self.Graph()
        if self.opts.write:
            self.WriteData()

    def RunAnalysis(self):
        r"""Run analysis
        """
        if self.verbose: print("GromacsAnalysis::RunAnalysis")
        # Check if we want to checkloadanalyze
        force_analyze = self.opts.force
        if not self.CheckLoadAnalyze(os.path.join(self.path, self.hd5_name),
                                     os.path.join(self.path, self.pickle_name),
                                     force_analyze):
            if self.verbose: print(f"  Loading previously analyzed data")
        else:
            if self.verbose: print(f"  Analyzing data")
            self.start_time = time.time()
            self.Analyze()
            print("--- %s seconds ---" % (time.time() - self.start_time))
        if self.verbose: print("GromacsAnalysis::RunAnalysis return")

    def CheckLoadAnalyze(self, file_path_pandas, file_path_pickle, force_analyze = False):
        r""" Check if the data can be loaded from HD5 and pickle files
        """
        if self.verbose: print("GromacsAnalysis::CheckLoadAnalyze")
        if self.verbose: print(f"  Forcing analysis: {force_analyze}")
        if os.path.isfile(file_path_pandas) and os.path.isfile(file_path_pickle) and not force_analyze:
            if self.verbose: print(f"  Found file(s), attempting load")
            # Skip analysis and load
            try:
                self.LoadData(file_path_pandas, file_path_pickle)
                if self.verbose: print("GromacsAnalysis::CheckLoadAnalyze return")
                return False
            except EOFError: return False
            except: raise
        else:
            if self.verbose: print("  Did not find file (or forcing load)")
            if self.verbose: print("GromacsAnalysis::CheckLoadAnalyze return")
            return True

    def Analyze(self):
        r""" Run analysis on a single simulation
        """
        if self.verbose: print("GromacsAnalysis::Analyze")

        # Create a reader depending on what we are reading in
        self.filename_structure = os.path.join(self.path, self.opts.structure)
        self.filename_trajectory = os.path.join(self.path, self.opts.trajectory)
        self.filename_gromacs = os.path.join(self.path, self.opts.gromacs)

        self.AnalyzeTrajectory()
        self.AnalyzeCurvature()
        self.AnalyzeDSSP()

        # Put the dataframes together
        self.master_time_df = pd.concat(self.time_dfs, axis=1)

        if self.verbose: print("GromacsAnalysis::Analyze return")

    def AnalyzeTrajectory(self):
        r""" Analyze trajectory information, unwrapped coordinates

        `AnalyzeTrajectory` is for all time-dependent analyses carried out in MDAnalysis
        """
        if self.verbose: print("GromacsAnalysis::AnalyzeTrajectory")
        # Create the universe
        traj_universe = mda.Universe(self.filename_structure, self.filename_trajectory)

        # Get the residues that are in the lipids
        resnames = set(np.unique(traj_universe.atoms.residues.resnames))
        common_lipids = resnames & self.lipid_set

        select_com_dict = {}
        lipid_selection = ''
        # Create selection logic for lipids
        for rname in common_lipids:
            select_com_dict[rname] = 'resname '+ rname
            lipid_selection = lipid_selection + select_com_dict[rname] + ' or '

        # Remove trailing or
        lipid_selection = lipid_selection[:-3]
        if self.verbose: print(f"  Lipid selection: {lipid_selection}")

        # Now get all the atoms that we want to analyze, also, get the COM structure set up
        atom_com_dict = {}
        lipid_com_dict = {}
        for key,val in select_com_dict.items():
            atom_com_dict[key] = traj_universe.select_atoms(val)
            lipid_com_dict[key] = []
        lipid_atoms =   traj_universe.select_atoms(lipid_selection)
        helix_atoms =   traj_universe.select_atoms('protein')
        not_protein =   traj_universe.select_atoms('not protein')
        not_solvent =   traj_universe.select_atoms(lipid_selection + ' or protein')
        solvent =       traj_universe.select_atoms('not (' + lipid_selection + ' or protein)')

        # Try to get the head groups of the two leaflets for analysis too
        # Has an implicit cutoff at 15.0
        L = leaf.LeafletFinder(traj_universe, 'name P*')
        leaflet0 = L.groups(0)
        leaflet1 = L.groups(1)

        # Unwrap/wrap the protein so that it doesn't have issues with the PBC
        transforms = [trans.unwrap(helix_atoms),
                      trans.unwrap(lipid_atoms),
                      trans.center_in_box(lipid_atoms, wrap=True),
                      trans.wrap(solvent),
                      trans.wrap(lipid_atoms)]
        traj_universe.trajectory.add_transformations(*transforms)

        # Do the analysis on the trajectory
        times = []
        lipid_com = []
        helix_com = []
        leaflet0_com = []
        leaflet1_com = []
        unit_cell = []
        # Wrap in a nice progress bar for ourselves, yay
        for ts in ProgressBar(traj_universe.trajectory):
            # Get the time of the snapshot
            times.append(traj_universe.trajectory.time)

            # Get the center of mass and store it
            for key,_ in lipid_com_dict.items():
                lipid_com_dict[key].append(atom_com_dict[key].center_of_mass())

            # Get the general center of mass as well
            lipid_com.append(lipid_atoms.center_of_mass())
            helix_com.append(helix_atoms.center_of_mass())
            leaflet0_com.append(leaflet0.center_of_mass())
            leaflet1_com.append(leaflet1.center_of_mass())

            # Get the unit cell as well
            unit_cell.append(traj_universe.dimensions)

        # Save off the times for other uses!
        self.times = times

        # Run a helix analysis
        self.helix_analysis = hel.HELANAL(traj_universe, select='name CA and resid 1-18').run()

        # Split up into a pandas dataframe for viewing
        # Get all of the lipid center of mass to store in the dataframe
        for key,val in lipid_com_dict.items():
            xname = key + '_x'
            yname = key + '_y'
            zname = key + '_z'
            df = pd.DataFrame(val, columns = [xname, yname, zname], index=times)
            df.index.name = 'Time(ps)'
            self.time_dfs.append(df)
        # Get the aggregated information on the lipids/leaflets/protein
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
        # Extract the helix parameters that we want from the helix_analysis
        global_tilt_df = pd.DataFrame(self.helix_analysis.results.global_tilts, columns = ['global_tilt'], index=times)
        global_tilt_df.index.name = 'Time(ps)'

        self.time_dfs.append(lipid_com_df)
        self.time_dfs.append(helix_com_df)
        self.time_dfs.append(leaflet0_com_df)
        self.time_dfs.append(leaflet1_com_df)
        self.time_dfs.append(unit_cell_df)
        self.time_dfs.append(global_tilt_df)
        
        if self.verbose: print("GromacsAnalysis::AnalyzeTrajectory return")

    def AnalyzeDSSP(self):
        r""" Analyze DSSP parameters for secondary structure
        """
        if self.verbose: print(f"GromacsAnalysis::AnalyzeDSSP")
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

        self.time_dfs.append(helicity_df)

        if self.verbose: print(f"GromacsAnalysis::AnalyzeDSSP return")

    def AnalyzeCurvature(self):
        r""" Analyze the curvature of the membrane
        """
        if self.verbose: print("GromacsAnalysis::AnalyzeCurvature")

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
        self.curvature_upper_leaflet = MembraneCurvature(universe,
                                                         select = upper_string,
                                                         n_x_bins = nbins,
                                                         n_y_bins = nbins,
                                                         x_range = x_dim,
                                                         y_range = y_dim,
                                                         wrap = True).run()
        self.curvature_lower_leaflet = MembraneCurvature(universe,
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

        # We have to go to a lot of trouble to save this in a pandas dataframe, but might as well, as then we can read it back in
        # in an easier setup. Basically, just iterate over the combinations of the columns, naming them as we go
        # Generate column names
        data_name       = ['z_surface', 'mean_curvature', 'gaussian_curvature']
        data_name_avg   = ['avg_z_surface', 'avg_mean_curvature', 'avg_gaussian_curvature']
        curvature_dict = {}
        for lf in leaflets:
            for dname in np.concatenate((data_name, data_name_avg), axis=0):
                for idx in np.arange(nbins):
                    for idy in np.arange(nbins):
                        full_name = lf + '_' + dname + '_' + str(idx) + '_' + str(idy)
                        curvature_dict[full_name] = []
        curvature_dict['curvature_helix_x'] = []
        curvature_dict['curvature_helix_y'] = []
        curvature_dict['curvature_helix_z'] = []

        # Generate the entries to the columns
        nframes = self.curvature_upper_leaflet.results.z_surface.shape[0]
        for idx in np.arange(nframes):
            for lf in leaflets:
                if lf == 'Upper':
                    mresults = self.curvature_upper_leaflet.results
                elif lf == 'Lower':
                    mresults = self.curvature_lower_leaflet.results

                for dname in data_name:
                    if dname == 'z_surface':
                        mresults_grid = mresults.z_surface
                    elif dname == 'mean_curvature':
                        mresults_grid = mresults.mean
                    elif dname == 'gaussian_curvature':
                        mresults_grid = mresults.gaussian

                    for xx in np.arange(nbins):
                        for yy in np.arange(nbins):
                            full_name = lf + '_' + dname + '_' + str(xx) + '_' + str(yy)
                            curvature_dict[full_name].append(mresults_grid[idx][xx][yy])

            curvature_dict['curvature_helix_x'].append(helix_com_new[idx][0])
            curvature_dict['curvature_helix_y'].append(helix_com_new[idx][1])
            curvature_dict['curvature_helix_z'].append(helix_com_new[idx][2])

        # Now just get the average values and store them too
        for lf in leaflets:
            if lf == 'Upper':
                mresults = self.curvature_upper_leaflet.results
            elif lf == 'Lower':
                mresults = self.curvature_lower_leaflet.results

            for dname in data_name_avg:
                if dname == 'avg_z_surface':
                    mresults_grid = mresults.average_z_surface
                elif dname == 'avg_mean_curvature':
                    mresults_grid = mresults.average_mean
                elif dname == 'avg_gaussian_curvature':
                    mresults_grid = mresults.average_gaussian

                for xx in np.arange(nbins):
                    for yy in np.arange(nbins):
                        full_name = lf + '_' + dname + '_' + str(xx) + '_' + str(yy)
                        curvature_dict[full_name].append(mresults_grid[xx][yy])

        # Save off the information as a dataframe
        curvature_df = pd.DataFrame(curvature_dict, index = self.times) 
        curvature_df.index.name = 'Time(ps)'
        self.time_dfs.append(curvature_df)

    def Graph(self):
        r""" Graph results
        """
        if self.verbose: print(f"GromacsAnalysis::Graph")

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
        fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
        self.graph_helicity(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_helicity.pdf', dpi=fig.dpi)

        # Plot helix information
        fig, axarr = plt.subplots(2, 1, figsize = (15, 10))
        self.graph_helix_analysis(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_helix_analysis.pdf', dpi=fig.dpi)

        # Plot the global tilts
        fig, axarr = plt.subplots(1, 1, figsize = (15, 10))
        self.graph_global_tilt(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_global_tilt.pdf', dpi=fig.dpi)

        # Plot the mean location of the upper leaflet
        fig, axarr = plt.subplots(1, 2, figsize = (15, 10))
        self.graph_avg_z_surface(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_avg_zsurf.pdf', dpi=fig.dpi)

        # Plot the average mean curvature
        fig, axarr = plt.subplots(1, 2, figsize = (15, 10))
        self.graph_avg_mean_curvature(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_avg_meancurv.pdf', dpi = fig.dpi)

        # Plot the average gaussian curvature
        fig, axarr = plt.subplots(1, 2, figsize = (15, 10))
        self.graph_avg_gaussian_curvature(axarr)
        fig.tight_layout()
        fig.savefig('gromacs_avg_gausscurv.pdf', dpi = fig.dpi)

        if self.verbose: print(f"GromacsAnalysis::Graph return")

    def WriteData(self):
        r""" Write the data into either CSV files for pandas or pickles for internal objects
        """
        if self.verbose: print(f"GromacsAnalysis::WriteData")
        # Dump the CSV files
        hd5_filename = os.path.join(self.path, self.name + '.h5')
        self.master_time_df.to_hdf(hd5_filename, key='master_time_df', mode='w')

        if self.verbose:
            print(self.master_time_df)

        # Dump the pickle files
        pkl_filename = os.path.join(self.path, self.name + '.pickle')
        with open(pkl_filename, 'wb') as f:
            pickle.dump(self.helix_analysis,f)

        if self.verbose: print(f"GromacsAnalysis::WriteData return")

    def LoadData(self, file_path_pandas, file_path_helixanalysis):
        r""" Load the data for pandas and pickles
        """
        if self.verbose: print(f"GromacsAnalysis::LoadData")

        self.master_time_df = pd.read_hdf(file_path_pandas)
        with open(file_path_helixanalysis, 'rb') as f:
            self.helix_analysis = pickle.load(f)

        if self.verbose: print(f"GromacsAnalysis::LoadData return")

    def graph_zdist_com(self, axarr):
        r""" Plot the Z distance between the lipid and helix COM
        """
        z = np.abs(self.master_time_df['helix_z'] - self.master_time_df['lipid_z'])
        axarr.plot(z)
        axarr.set_xlabel('Frame')
        axarr.set_ylabel('Z distance (Angstroms)')

    def graph_zpos_com(self, axarr):
        r""" Plot the Z positino of the lipid head groups and the protein
        """
        z_protein = self.master_time_df['helix_z']
        z_leaf0 = self.master_time_df['leaflet0_z']
        z_leaf1 = self.master_time_df['leaflet1_z']
        z_lipid = self.master_time_df['lipid_z']
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
        helicity = self.master_time_df['helicity']
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

    def graph_global_tilt(self, axarr):
        r""" Plot the global tilt angle
        """
        global_tilt = self.master_time_df['global_tilt']
        axarr.plot(global_tilt, color = 'k')
        axarr.set_xlabel('Frame')
        axarr.set_ylabel('Global Tilt (deg)')

    def graph_avg_z_surface(self, axarr):
        r""" Plot the average z surface
        """
        # Unpack the curvature dataframe
        lower_average_z_surface = unpack_curvature(self.master_time_df, 'Lower', 'avg_z_surface')
        upper_average_z_surface = unpack_curvature(self.master_time_df, 'Upper', 'avg_z_surface')

        surfaces = [lower_average_z_surface[0][:][:],
                    upper_average_z_surface[0][:][:]]
        leaflets = ['Lower', 'Upper']
        for ax, surfs, lf in zip(axarr, surfaces, leaflets):
            im = ax.imshow(surfs, interpolation = 'gaussian', cmap = 'YlGnBu', origin = 'lower')
            ax.set_aspect('equal')
            ax.set_title('{} Leaflet'.format(lf))
            cbar = plt.colorbar(im, ticks = [surfs.min(), surfs.max()], orientation='horizontal', ax=ax, shrink=0.7, aspect=10, pad=0.05)
            cbar.set_ticklabels([int(surfs.min()), int(surfs.max())])
            cbar.ax.tick_params(width=0.5)
            cbar.set_label("Height lipid headgroups (${\AA}$)", labelpad=2)

    def graph_avg_mean_curvature(self, axarr):
        r""" Plot the average mean curvature
        """
        # Unpack the curvature dataframe
        lower_average_mean_curvature = unpack_curvature(self.master_time_df, 'Lower', 'avg_mean_curvature')
        upper_average_mean_curvature = unpack_curvature(self.master_time_df, 'Upper', 'avg_mean_curvature')

        surfaces = [lower_average_mean_curvature[0][:][:],
                    upper_average_mean_curvature[0][:][:]]
        leaflets = ['Lower', 'Upper']
        for ax, surfs, lf in zip(axarr, surfaces, leaflets):
            im = ax.imshow(surfs, interpolation = 'gaussian', cmap = 'bwr', origin = 'lower')
            ax.set_aspect('equal')
            ax.set_title('{} Leaflet'.format(lf))
            cbar = plt.colorbar(im, ticks = [surfs.min(), surfs.max()], orientation='horizontal', ax=ax, shrink=0.7, aspect=10, pad=0.05)
            cbar.set_ticklabels([int(surfs.min()), int(surfs.max())])
            cbar.ax.tick_params(width=0.5)
            cbar.set_label("$H$ (${\AA}^{-1}$)", labelpad=2)

    def graph_avg_gaussian_curvature(self, axarr):
        r""" Plot the average mean curvature
        """
        # Unpack the curvature dataframe
        lower_average_gaussian_curvature = unpack_curvature(self.master_time_df, 'Lower', 'avg_gaussian_curvature')
        upper_average_gaussian_curvature = unpack_curvature(self.master_time_df, 'Upper', 'avg_gaussian_curvature')

        surfaces = [lower_average_gaussian_curvature[0][:][:],
                    upper_average_gaussian_curvature[0][:][:]]
        leaflets = ['Lower', 'Upper']
        for ax, surfs, lf in zip(axarr, surfaces, leaflets):
            im = ax.imshow(surfs, interpolation = 'gaussian', cmap = 'PiYG', origin = 'lower')
            ax.set_aspect('equal')
            ax.set_title('{} Leaflet'.format(lf))
            cbar = plt.colorbar(im, ticks = [surfs.min(), surfs.max()], orientation='horizontal', ax=ax, shrink=0.7, aspect=10, pad=0.05)
            cbar.set_ticklabels([int(surfs.min()), int(surfs.max())])
            cbar.ax.tick_params(width=0.5)
            cbar.set_label("$K$ (${\AA}^{-2}$)", labelpad=2)





##########################################
if __name__ == "__main__":
    opts = parse_args()
    x = GromacsAnalysis(opts)
