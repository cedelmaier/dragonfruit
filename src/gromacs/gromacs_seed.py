#!/usr/bin/env python3

# XXX: Put a license here

"""Main analysis script for membranes with AH domains"""

import logging
import pickle
import os
import sys
import time
import yaml

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

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from stylelib.common_styles import septin_runs_stl
from seed_base import SeedBase
from seed_graph_funcs import *

class GromacsSeed(SeedBase):
    def __init__(self, path, opts):
        # Try to set up the logging facilities correctly
        mda_logger = logging.getLogger('MDAnalysis')
        if mda_logger.hasHandlers():
            mda_logger.handlers = []
        mda.start_logging()
        self.logger = logging.getLogger('MDAnalysis')

        self.verbose = opts.verbose
        if self.verbose: print("GromacsSeed::__init__")

        SeedBase.__init__(self, path, opts)

        self.ReadData()

        # Set up the named versions of files at the end
        self.hd5_name       = self.name + '.h5'
        self.pickle_name    = self.name + '.pickle'
        self.time_dfs = []
        self.lipid_set = {"DOPC",
                          "PLPI",
                          "CHL1",
                          "POPC",
                          "DOPE",
                          "DOPS",
                          "SAPI24",
                          "SAPI25"}

        if self.verbose: print("GromacsSeed::__init__ return")

    def ReadData(self):
        r""" Read the configuration from the YAML file for analysis
        """
        if self.verbose: print("GromacsSeed::ReadData")

        # Configuration
        self.structure_file     = self.default_yaml['structure']
        self.trajectory_file    = self.default_yaml['trajectory']
        self.gromacs_file       = self.default_yaml['gromacs']

        if self.verbose: print("GromacsSeed::ReadData return")

    def PrintInformation(self):
        r""" Print parameter information for GromacsSeed
        """
        print(f"--------")
        print(f"Seed {self.name}")
        print(f"--------")
        print(f"File information")
        print(f"Structure file      = {self.structure_file}")
        print(f"Trajectory file     = {self.trajectory_file}")
        print(f"Gromacs file        = {self.gromacs_file}")
        print(f"--------")

    def CheckLoadAnalyze(self, file_path_pandas, file_path_pickle, force_analyze = False):
        r""" Check if the data can be loaded from HD5 and pickle files
        """
        if self.verbose: print("GromacsSeed::CheckLoadAnalyze")
        if self.verbose: print(f"  Forcing analysis: {force_analyze}")
        if os.path.isfile(file_path_pandas) and os.path.isfile(file_path_pickle) and not force_analyze:
            if self.verbose: print(f"  Found file(s), attempting load")
            # Skip analysis and load
            try:
                self.LoadData(file_path_pandas, file_path_pickle)
                if self.verbose: print("GromacsSeed::CheckLoadAnalyze return")
                return False
            except EOFError: return False
            except: raise
        else:
            if self.verbose: print("  Did not find file (or forcing load)")
            if self.verbose: print("GromacsSeed::CheckLoadAnalyze return")
            return True

    def LoadData(self, file_path_pandas, file_path_helixanalysis):
        r""" Load the data from pandas and pickle
        """
        if self.verbose: print(f"GromacsSeed::LoadData")

        self.master_time_df = pd.read_hdf(file_path_pandas)
        with open(file_path_helixanalysis, 'rb') as f:
            self.helix_analysis = pickle.load(f)

        if self.verbose: print(f"GromacsSeed::LoadData return")

    def Analyze(self, force_analyze = False):
        r""" Analysis of a single gromacs simulation seed
        """
        if self.verbose: print("GromacsSeed::Analyze")
        # Check the force flag
        if self.opts.force:
            force_analyze = True

        # Check if we are loading or analyzing the information
        if not self.CheckLoadAnalyze(os.path.join(self.path, self.hd5_name),
                                     os.path.join(self.path, self.pickle_name),
                                     force_analyze):
            if self.verbose: print("  Loading previously analyzed data")
            if self.verbose: print("GromacsSeed::Analyze return")
            return False

        if self.verbose: print("  Analyzing data")
        
        self.PrintInformation()

        self.start_time = time.time()

        # Main analysis loop, or subloops
        # Create a reader depending on what we are reading in
        self.filename_structure     = os.path.join(self.path, self.structure_file)
        self.filename_trajectory    = os.path.join(self.path, self.trajectory_file)
        self.filename_gromacs       = os.path.join(self.path, self.gromacs_file)

        self.AnalyzeTrajectory()
        self.AnalyzeCurvature()
        self.AnalyzeDSSP()

        # Put the master time dataframe together
        self.master_time_df = pd.concat(self.time_dfs, axis = 1)

        # Save the data
        self.WriteData()

        if self.verbose: print("---- %s seconds ----" % (time.time() - self.start_time))
        if self.verbose: print("GromacsSeed::Analyze return")

    def AnalyzeTrajectory(self):
        r""" Analyze trajectory information, unwrapped coordinates

        `AnalyzeTrajectory` is for all time-dependent analyses carried out in MDAnalysis
        """
        if self.verbose: print("GromacsSeed::AnalyzeTrajectory")
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
        p_dipole_list = []
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

            # Calculate the dipole moment
            nhelix = len(helix_atoms.positions[:,2])
            p_r_dipole = np.zeros((nhelix, 3))
            p_dipole = np.array([0.0, 0.0, 0.0])
            for iatom in range(nhelix):
                atom_name   = helix_atoms.atoms.names[iatom]
                atom_charge = helix_atoms.atoms.charges[iatom]
                relative_position = (helix_atoms.positions[iatom,:] - helix_atoms.center_of_mass())
                p_r_dipole[iatom] = relative_position * atom_charge
                p_dipole += p_r_dipole[iatom]
            p_dipole_list.append(p_dipole)

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
        global_tilt_df = pd.DataFrame(self.helix_analysis.results.global_tilts, columns = ['global_tilt'], index = times)
        global_tilt_df.index.name = 'Time(ps)'
        global_axis_df = pd.DataFrame(self.helix_analysis.results.global_axis, columns = ['helix_global_axis_x', 'helix_global_axis_y', 'helix_global_axis_z'], index = times)
        global_axis_df.index.name = 'Time(ps)'
        # What about the dipole moment of the helix?
        p_dipole_df = pd.DataFrame(p_dipole_list, columns = ['p_dipole_x', 'p_dipole_y', 'p_dipole_z'], index = times)
        p_dipole_df.index.name = 'Time(ps)'

        self.time_dfs.append(lipid_com_df)
        self.time_dfs.append(helix_com_df)
        self.time_dfs.append(leaflet0_com_df)
        self.time_dfs.append(leaflet1_com_df)
        self.time_dfs.append(unit_cell_df)
        self.time_dfs.append(global_tilt_df)
        self.time_dfs.append(global_axis_df)
        self.time_dfs.append(p_dipole_df)
        
        if self.verbose: print("GromacsSeed::AnalyzeTrajectory return")

    def AnalyzeDSSP(self):
        r""" Analyze DSSP parameters for secondary structure
        """
        if self.verbose: print(f"GromacsSeed::AnalyzeDSSP")
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

        if self.verbose: print(f"GromacsSeed::AnalyzeDSSP return")

    def AnalyzeCurvature(self):
        r""" Analyze the curvature of the membrane
        """
        if self.verbose: print("GromacsSeed::AnalyzeCurvature")

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

        if self.verbose: print("GromacsSeed::AnalyzeCurvature return")

    def WriteData(self):
        r""" Write the data to HD5 and pickle files
        """
        if self.verbose: print(f"GromacsSeed::WriteData")

        # Dump the HD5 files
        hd5_filename = os.path.join(self.path, self.hd5_name)
        self.master_time_df.to_hdf(hd5_filename, key='master_time_df', mode='w')

        if self.verbose:
            print(self.master_time_df)

        # Dump the pickle file(s)
        pickle_filename = os.path.join(self.path, self.pickle_name)
        with open(pickle_filename, 'wb') as f:
            pickle.dump(self.helix_analysis,f)

        if self.verbose: print(f"GromacsSeed::WriteData return")

    def GraphZdist(self, ax, color = 'b'):
        r""" Plot the absolute Z position of the helix
        """
        graph_seed_zdist(self, ax, color = color)

    def GraphZpos(self, ax, color = 'b'):
        r""" Plot the Z position of the lipid headgroups and protein
        """
        graph_seed_zpos_wheads(self, ax, color = color)

    def GraphHelicity(self, ax, color = 'b'):
        r""" Plot the fractional helicity of the helix
        """
        graph_seed_helicity(self, ax, color = 'b') 

    def GraphGlobalTilt(self, ax, color = 'b'):
        r""" Plot the global tilt of the helix
        """
        graph_seed_globaltilt(self, ax, color = 'b')

    def GraphPDipoleTilt(self, ax, color = 'b'):
        r""" Plot the tilt of the helix dipole with respect to the Z axis
        """
        graph_seed_pdipoletilt(self, ax, color = 'b')

    def GraphHelixPDipoleAngle(self, ax, color = 'b'):
        r""" Plot the angle between the helix vector and the p_dipole vector
        """
        graph_seed_helixpdipoleangle(self, ax, color = 'b')

    def GraphPDipoleMoment(self, ax, color = 'b'):
        r""" Plot the angle between the helix vector and the p_dipole vector
        """
        graph_seed_pdipolemoment(self, ax, color = 'b')

#    def graph_helix_analysis(self, axarr):
#        r""" Plot results of the helix analysis
#        """
#        axarr[0].plot(self.helix_analysis.results.local_twists.mean(axis=1))
#        axarr[0].set_xlabel('Frame')
#        axarr[0].set_ylabel('Average twist (degrees)')
#
#        axarr[1].plot(self.helix_analysis.results.local_nres_per_turn.mean(axis=1))
#        axarr[1].set_xlabel('Frame')
#        axarr[1].set_ylabel('Average residues per turn')

    def graph_avg_z_surface(self, axarr):
        r""" Plot the average z surface
        """
        # Unpack the curvature dataframe
        lower_average_z_surface = self.unpack_curvature(self.master_time_df, 'Lower', 'avg_z_surface')
        upper_average_z_surface = self.unpack_curvature(self.master_time_df, 'Upper', 'avg_z_surface')

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
        lower_average_mean_curvature = self.unpack_curvature(self.master_time_df, 'Lower', 'avg_mean_curvature')
        upper_average_mean_curvature = self.unpack_curvature(self.master_time_df, 'Upper', 'avg_mean_curvature')

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
        lower_average_gaussian_curvature = self.unpack_curvature(self.master_time_df, 'Lower', 'avg_gaussian_curvature')
        upper_average_gaussian_curvature = self.unpack_curvature(self.master_time_df, 'Upper', 'avg_gaussian_curvature')

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

    def unpack_curvature(self, df, leaflet, measurement):
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

