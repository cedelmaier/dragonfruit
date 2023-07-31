#!/usr/bin/env python3

# XXX: Put a license here

"""Main analysis script for membranes with AH domains"""

import logging
import pickle
import os
import re
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
from dragonfruit_io import create_chargeless_topology
from gromacs_electrostatics import GromacsElectrostaticsAnalyzer

class GromacsSeed(SeedBase):
    def __init__(self, path, opts):
        # Try to set up the logging facilities correctly
        mda_logger = logging.getLogger('MDAnalysis')
        if mda_logger.hasHandlers():
            mda_logger.handlers = []
        mda.start_logging()
        self.logger = logging.getLogger('MDAnalysis')

        self.verbose = opts.verbose
        self.trace = opts.trace
        if self.verbose: print("GromacsSeed::__init__")

        SeedBase.__init__(self, path, opts)

        self.ReadData()

        # Set up the named versions of files at the end
        self.hd5_name               = self.name + '.h5'
        self.electrostatics_name    = self.name + '_electrostatics.h5'
        self.forces_name            = self.name + '_forces.h5'
        self.pickle_name            = self.name + '.pickle'
        self.time_dfs = []
        self.new_electrostatics_dfs = []
        self.lipid_set = {"DOPC",
                          "PLPI",
                          "CHL1",
                          "POPC",
                          "DOPE",
                          "DOPS",
                          "SAPI24",
                          "SAPI25"}

        # Create any subclasses that we might need
        self.gromacs_electrostatics = GromacsElectrostaticsAnalyzer(path, opts, self.structure_file)

        if self.verbose: print("GromacsSeed::__init__ return")

    def ReadData(self):
        r""" Read the configuration from the YAML file for analysis
        """
        if self.verbose: print("GromacsSeed::ReadData")

        # Configuration
        self.structure_file     = self.default_yaml['structure']
        self.trajectory_file    = self.default_yaml['trajectory']
        self.gromacs_file       = self.default_yaml['gromacs']

        # Do some special stuff if we detect a plumed section, this is to not break backwards compatibility
        if 'plumed' in self.default_yaml:
            self.plumed_ref     = self.default_yaml['plumed']['reference_file']
            self.plumed_file    = self.default_yaml['plumed']['plumed_file']
            self.alpha_residues = self.default_yaml['plumed']['alpha_residues']
            self.lipid_P_atoms  = self.default_yaml['plumed']['lipid_P_atoms']
            self.helix_C_atoms  = self.default_yaml['plumed']['helix_C_atoms']
            self.colvar_file    = self.default_yaml['plumed']['colvar_file']
            self.mcfile         = self.default_yaml['plumed']['masscharge_file']

            # Check for the specializations of the whole helix
            if 'wholehelix_C_atoms' in self.default_yaml['plumed']:
                self.wholehelix_C_atoms     = self.default_yaml['plumed']['wholehelix_C_atoms']
                self.wholehelix_residues    = self.default_yaml['plumed']['wholehelix_residues']
            else:
                self.wholehelix_C_atoms     = self.helix_C_atoms
                self.wholehelix_residues    = self.alpha_residues


        # Do something if we see information for the N and C terminals
        if 'N_term_Calpha' in self.default_yaml:
            self.N_term_Calpha  = self.default_yaml['N_term_Calpha']
            self.C_term_Calpha  = self.default_yaml['C_term_Calpha']

        # Check to see if we are doing electrostatics or not
        self.do_electrostatics = True
        if 'do_electrostatics' in self.default_yaml:
            self.do_electrostatics = self.default_yaml['do_electrostatics']

        # Harcoded information
        self.zbin_size          = 0.5   # Size of z-bins for calculations
        self.zbin_range         = 40.0  # Range of Z bins for calculations

        # We keep track of the histogram information ourself due to various reasons
        self.zbin_nbins = np.int32(self.zbin_range / self.zbin_size)
        self.zbin_edges = np.linspace(-self.zbin_range, self.zbin_range, self.zbin_nbins+1)
        self.zbin_mids = moving_average(self.zbin_edges)

        if self.verbose: print("GromacsSeed::ReadData return")

    def PrintInformation(self):
        r""" Print parameter information for GromacsSeed
        """
        print(f"--------")
        print(f"Seed {self.name}")
        print(f"--------")
        print(f"File information")
        print(f"Path                = {self.path}")
        print(f"Structure file      = {self.structure_file}")
        print(f"Trajectory file     = {self.trajectory_file}")
        print(f"Gromacs file        = {self.gromacs_file}")
        print(f"--------")
        print(f"Electrostatic information")
        print(f"Do electrostatics   = {self.do_electrostatics}")
        print(f"Z-bin size          = {self.zbin_size}")
        print(f"Z-range             = {self.zbin_range}")
        print(f"Nbins EM            = {self.zbin_nbins}")
        if 'plumed' in self.default_yaml:
            print(f"--------")
            print(f"Plumed information")
            print(f"PLUMED reference    = {self.plumed_ref}")
            print(f"PLUMED file         = {self.plumed_file}")
            print(f"Colvar file         = {self.colvar_file}")
            print(f"Mass/charge file    = {self.mcfile}")
            print(f"Alpha residues      = {self.alpha_residues}")
            print(f"Lipid P atoms       = {self.lipid_P_atoms}")
            print(f"Helix C atoms       = {self.helix_C_atoms}")
            print(f"Whole helix C atoms = {self.wholehelix_C_atoms}")
            print(f"Whole helix residue = {self.wholehelix_residues}")
        if 'N_term_Calpha' in self.default_yaml:
            print(f"--------")
            print(f"N and C term information")
            print(f"N term Calphas      = {self.N_term_Calpha}")
            print(f"C term Calphas      = {self.C_term_Calpha}")

        print(f"--------")

    def CheckLoadAnalyze(self, force_analyze = False):
        r""" Check if the data can be loaded from HD5 and pickle files
        """
        if self.verbose: print("GromacsSeed::CheckLoadAnalyze")
        if self.verbose: print(f"  Forcing analysis: {force_analyze}")

        # Set up the file paths for ourselves
        file_path_pandas = os.path.join(self.path, self.hd5_name)
        file_path_electrostatics = os.path.join(self.path, self.electrostatics_name)
        file_path_forces = os.path.join(self.path, self.forces_name)
        file_path_pickle = os.path.join(self.path, self.pickle_name)

        # If there was a change in name at some point, try to account for that, as if we cannot find the .h5 files
        # we might have a different name
        seed_name_regex = re.compile('N[0-9]*.h5')
        current_h5_name = self.name + ".h5"
        for root, dirs, files in os.walk(self.path):
            for file in files:
                if seed_name_regex.match(file):
                    if file != current_h5_name:
                        print(f"  Found a different .h5 file, setting new names")
                        filebase = file.split('.')[0]
                        file_path_pandas = os.path.join(self.path, filebase + ".h5")
                        file_path_electrostatics = os.path.join(self.path, filebase + "_electrostatics.h5")
                        file_path_forces = os.path.join(self.path, filebase + "_forces.h5")
                        file_path_pickle = os.path.join(self.path, filebase + ".pickle")
                        break

        if os.path.isfile(file_path_pandas) and os.path.isfile(file_path_pickle) and not force_analyze:
            if self.verbose: print(f"  Found file(s), attempting load")
            # Skip analysis and load
            try:
                self.LoadData(file_path_pandas, file_path_electrostatics, file_path_forces, file_path_pickle)
                if self.verbose: print("GromacsSeed::CheckLoadAnalyze return")
                return False
            except EOFError: return False
            except: raise
        else:
            if self.verbose: print("  Did not find file (or forcing load)")
            if self.verbose: print("GromacsSeed::CheckLoadAnalyze return")
            return True

    def LoadData(self, file_path_pandas, file_path_electrostatics, file_path_forces, file_path_helixanalysis):
        r""" Load the data from pandas and pickle
        """
        if self.verbose: print(f"GromacsSeed::LoadData")

        self.master_time_df = pd.read_hdf(file_path_pandas)
        with open(file_path_helixanalysis, 'rb') as f:
            self.helix_analysis = pickle.load(f)

        if self.do_electrostatics:
            self.gromacs_electrostatics.electrostatic_df = pd.read_hdf(file_path_electrostatics)
            self.master_forces_df = pd.read_hdf(file_path_forces)

        if self.verbose: print(f"GromacsSeed::LoadData return")

    def Analyze(self, force_analyze = False):
        r""" Analysis of a single gromacs simulation seed
        """
        if self.verbose: print("GromacsSeed::Analyze")
        print(f"Seed {self.name} reporting for duty!")
        # Check the force flag
        if self.opts.force:
            force_analyze = True

        # Check if we are loading or analyzing the information
        if not self.CheckLoadAnalyze(force_analyze):
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

        self.PrepareAnalysis()

        # Order here does matter, as need electrostatics before trajectory
        if self.do_electrostatics:
            self.AnalyzeElectrostatics()
        self.AnalyzeTrajectory()
        self.AnalyzeCurvature()
        self.AnalyzeDSSP()
        self.AnalyzeAlphaRMSD()

        # Put the master time dataframe together
        self.master_time_df = pd.concat(self.time_dfs, axis = 1)

        # Put the force calculations together
        if self.do_electrostatics:
            self.master_forces_df = pd.concat(self.new_electrostatics_dfs, axis = 1)

        # Save the data
        self.WriteData()

        if self.verbose: print("---- %s seconds ----" % (time.time() - self.start_time))
        if self.verbose: print("GromacsSeed::Analyze return")

    def PrepareAnalysis(self):
        r""" Prepare for the analysis by pre-calculating certain quantities
        """
        if self.verbose: print("GromacsSeed::PrepareAnalysis")

        # Create the universe
        traj_universe = mda.Universe(self.filename_structure, self.filename_trajectory)

        # Get the residues that are in the lipids
        resnames = set(np.unique(traj_universe.atoms.residues.resnames))
        self.common_lipids = resnames & self.lipid_set

        self.select_com_dict = {}
        self.lipid_selection = ''
        # Create selection logic for lipids
        for rname in self.common_lipids:
            self.select_com_dict[rname] = 'resname '+ rname
            self.lipid_selection = self.lipid_selection + self.select_com_dict[rname] + ' or '

        # Remove trailing or
        self.lipid_selection = self.lipid_selection[:-3]
        if self.verbose: print(f"  Lipid selection: {self.lipid_selection}")

        if self.verbose: print("GromacsSeed::PrepareAnalysis return")

    def AnalyzeTrajectory(self):
        r""" Analyze trajectory information, unwrapped coordinates

        `AnalyzeTrajectory` is for all time-dependent analyses carried out in MDAnalysis
        """
        if self.verbose: print("GromacsSeed::AnalyzeTrajectory")

        # Creat the universe
        traj_universe = mda.Universe(self.filename_structure, self.filename_trajectory)

        # Now get all the atoms that we want to analyze, also, get the COM structure set up
        atom_com_dict = {}
        lipid_com_dict = {}
        for key,val in self.select_com_dict.items():
            atom_com_dict[key] = traj_universe.select_atoms(val)
            lipid_com_dict[key] = []
        lipid_atoms =   traj_universe.select_atoms(self.lipid_selection)
        helix_atoms =   traj_universe.select_atoms('protein')
        not_protein =   traj_universe.select_atoms('not protein')
        not_solvent =   traj_universe.select_atoms(self.lipid_selection + ' or protein')
        solvent =       traj_universe.select_atoms('not (' + self.lipid_selection + ' or protein)')

        n_term_atoms = None
        c_term_atoms = None
        if 'N_term_Calpha' in self.default_yaml:
            atom_indices = self.N_term_Calpha
            atom_indices_shifted = [i-1 for i in atom_indices]
            selection_string = 'id ' + ' '.join(map(str, atom_indices_shifted))
            n_term_atoms = traj_universe.select_atoms(selection_string)

            atom_indices = self.C_term_Calpha
            atom_indices_shifted = [i-1 for i in atom_indices]
            selection_string = 'id ' + ' '.join(map(str, atom_indices_shifted))
            c_term_atoms = traj_universe.select_atoms(selection_string)

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

        n_term_com = []
        c_term_com = []

        # Keep track of the running histograms
        self.histograms = {}
        self.histograms['System'] = np.zeros(len(self.zbin_mids))

        # Forces
        force_breakdown_types = ["all", "noq", "q"]
        force_calc_types = ["force", "moment", "torque"]
        forces_com = {} # Master dictionary for all of the COM forces we are going to look at here, yay us
        # Initialize the forces
        for force_calc in force_calc_types:
            forces_com[force_calc] = {}
            for force_type in force_breakdown_types:
                forces_com[force_calc][force_type] = []
        # Create a copy of the electorstatic_df
        if self.do_electrostatics:
            elec_df = self.gromacs_electrostatics.electrostatic_df.copy(deep = True)

        # Run a helix analysis first, as we need its information for figuring out the location
        # of charges in the reference frame of the helix
        self.helix_analysis = hel.HELANAL(traj_universe, select='name CA and resid 1-18').run()

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

            # Check to see if we have N and C term information
            if n_term_atoms:
                n_term_com.append(n_term_atoms.center_of_mass())
                c_term_com.append(c_term_atoms.center_of_mass())

            # Get the unit cell as well
            unit_cell.append(traj_universe.dimensions)

            # Calculate the dipole moment
            if self.do_electrostatics:
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

                ## Calculte the charge density of various quantitites in the system
                #print(f"Trying to figure out charge densities")
                ## Loop through the edges of the zbins and find the charge density in each
                #for ibin in np.arange(len(self.zbin_mids)):
                #    print(f"Edge: {self.zbin_edges[ibin]}, {self.zbin_edges[ibin+1]}")
                #sys.exit(1)

                # Calculate the center of mass force
                # Figure out if we are synchronized for the two dataframes, as one is resampled
                if traj_universe.trajectory.time in elec_df.index:
                    # Loop over the force type we are intersted in to extract columns
                    for force_type in force_breakdown_types:
                        elec_df_type = elec_df.filter(regex = f".*_{force_type}")

                        # Extract just the row we are interested in
                        elec_df_row = elec_df_type.loc[traj_universe.trajectory.time]
                        # Cast this as a reshaped array with the proper [atom, dimension] setup
                        atom_forces = elec_df_row.to_numpy(dtype=np.float64).reshape(nhelix, 3)

                        # Have to calculate the net force first, as ths is used in the torque calculations
                        f_net = np.zeros(3, dtype=np.float64)
                        for iatom in range(nhelix):
                            #atom_name = helix_atoms.atoms.names[iatom]
                            #r_i = (helix_atoms.positions[iatom,:] - helix_atoms.center_of_mass())
                            f_i = atom_forces[iatom,:]
                            f_net += f_i
                        # Reprocess and do the torque and moment
                        f_moment = np.zeros((3, 3), dtype=np.float64)
                        torque = np.zeros(3, dtype=np.float64)
                        for iatom in range(nhelix):
                            r_i = (helix_atoms.positions[iatom,:] - helix_atoms.center_of_mass())
                            f_i = atom_forces[iatom,:] - f_net
                            torque += np.cross(r_i, f_i)
                            f_moment += np.outer(r_i, f_i)

                        # Record the products in their places
                        forces_com["force"][force_type].append(f_net)
                        forces_com["torque"][force_type].append(torque)
                        forces_com["moment"][force_type].append(f_moment)

        # END of progress bar loop

        # Add the forces to the electrostatics dataframe
        if self.do_electrostatics:
            for force_type in force_breakdown_types:
                # Force and torque follow the same convention
                force_column_names = ["force_com_" + force_type + "_" + xval for xval in ["x", "y", "z"]]
                force_df = pd.DataFrame(forces_com["force"][force_type], columns = force_column_names, index = elec_df.index)
                torque_column_names = ["torque_com_" + force_type + "_" + xval for xval in ["x", "y", "z"]]
                torque_df = pd.DataFrame(forces_com["torque"][force_type], columns = torque_column_names, index = elec_df.index)

                # Moment is slightly harder, but still fine
                n_moment_frames = len(force_df.index)
                moment_dict = {}
                x_names = ["x", "y", "z"]
                for idx in np.arange(n_moment_frames):
                    for xx in np.arange(3):
                        for yy in np.arange(3):
                            full_moment_name = "moment_com_" + force_type + "_" + x_names[xx] + x_names[yy]
                            if full_moment_name not in moment_dict:
                                moment_dict[full_moment_name] = []
                            moment_dict[full_moment_name].append(forces_com["moment"][force_type][idx][xx][yy])
                moment_df = pd.DataFrame(moment_dict, index = elec_df.index)

                self.new_electrostatics_dfs.append(force_df)
                self.new_electrostatics_dfs.append(torque_df)
                self.new_electrostatics_dfs.append(moment_df)

        # Save off the times for other uses!
        self.times = times

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
        if self.do_electrostatics:
            p_dipole_df = pd.DataFrame(p_dipole_list, columns = ['p_dipole_x', 'p_dipole_y', 'p_dipole_z'], index = times)
            p_dipole_df.index.name = 'Time(ps)'
        # What if we have n and c terms?
        if n_term_atoms:
            n_term_df = pd.DataFrame(n_term_com, columns = ['nterm_com_x', 'nterm_com_y', 'nterm_com_z'], index = times)
            n_term_df.index.name = 'Time(ps)'
            c_term_df = pd.DataFrame(c_term_com, columns = ['cterm_com_x', 'cterm_com_y', 'cterm_com_z'], index = times)
            c_term_df.index.name = 'Time(ps)'

        self.time_dfs.append(lipid_com_df)
        self.time_dfs.append(helix_com_df)
        self.time_dfs.append(leaflet0_com_df)
        self.time_dfs.append(leaflet1_com_df)
        self.time_dfs.append(unit_cell_df)
        self.time_dfs.append(global_tilt_df)
        self.time_dfs.append(global_axis_df)
        if self.do_electrostatics:
            self.time_dfs.append(p_dipole_df)
        if n_term_atoms:
            self.time_dfs.append(n_term_df)
            self.time_dfs.append(c_term_df)
        
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

        self.helicity   = np.zeros(helicity_results.shape[0])
        self.hel_alpha  = np.zeros(helicity_results.shape[0])
        for iframe,val in enumerate(helicity_results):
            # Get the helicity of the combined states
            h_count = (val == 'H').sum() + (val == 'G').sum() + (val == 'I').sum()
            h_total = val.shape

            current_helicity = h_count / h_total
            self.helicity[iframe] = current_helicity[0]
            self.hel_alpha[iframe] = h_count

        # Add the helicity to the master DF
        helicity_df = pd.DataFrame(self.helicity, columns = ['helicity'], index=self.times)
        helicity_df.index.name = 'Time(ps)'
        hel_alpha_df = pd.DataFrame(self.hel_alpha, columns = ['alpha_dssp'], index = self.times)
        hel_alpha_df.index.name = 'Time(ps)'

        self.time_dfs.append(helicity_df)
        self.time_dfs.append(hel_alpha_df)

        if self.verbose: print(f"GromacsSeed::AnalyzeDSSP return")

    def AnalyzeAlphaRMSD(self):
        r""" Analyze alpha_rmsd measurement from PLUMED
        """
        if self.verbose: print(f"GromacsSeed::AnalyzeAlphaRMSD")

        # Have to create a plumed file and run it, seek inspiration from gromacs_electrostatics package
        # Write the plumed file
        self.WritePLUMEDAlphaRMSD()

        # Run the plumed file
        import subprocess
        plumed_p = subprocess.Popen(["plumed", "driver", "--mf_xtc", self.trajectory_file, "--plumed", self.plumed_file, "--kt", "2.494339", "--mc", self.mcfile], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        plumed_output, plumed_errors = plumed_p.communicate()
        plumed_p.wait()
        if self.trace: print(plumed_errors.decode("utf-8"))
        if self.trace: print(plumed_output.decode("utf-8"))

        # Now that we have a colvar file, read it in via plumed
        import plumed
        pdata = plumed.read_as_pandas("./{}".format(self.colvar_file))
        # Set the times to the time_dataframe stuff
        pdata.rename(columns={"dz": "plumed_dz", "alpha": "plumed_alpha"}, inplace=True)
        pdata['time'] = self.times
        pdata.set_index('time', inplace=True)
        pdata.index.name = 'Time(ps)'

        # Add to the dataframe list
        self.time_dfs.append(pdata)

        if self.verbose: print(f"GromacsSeed::AnalyzeAlphaRMSD return")

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

    def AnalyzeElectrostatics(self):
        r""" Analyze electrostatics in the simluation

        NOTE: This is complicated, and requires access to a gromacs executable, and may
        cause issues when run on supercomputers, for instance. Be careful with this, as it
        also might cause large amounts of data to be written if not cleaned up properly!
        """
        if self.verbose: print("GromacsSeed::AnalyzeElectrostatics")

        self.gromacs_electrostatics.CreateChargelessTopology(self.common_lipids)

        self.gromacs_electrostatics.GenerateForces("all", self.trajectory_file)
        self.gromacs_electrostatics.GenerateForces("noq", self.trajectory_file)

        self.gromacs_electrostatics.CalculateElectrostaticForces()

        if self.verbose: print("GromacsSeed::AnalyzeElectrostatics return")

    def WriteData(self):
        r""" Write the data to HD5 and pickle files
        """
        if self.verbose: print(f"GromacsSeed::WriteData")

        # Dump the HD5 files
        hd5_filename = os.path.join(self.path, self.hd5_name)
        self.master_time_df.to_hdf(hd5_filename, key='master_time_df', mode='w')
        if self.verbose: print(self.master_time_df)

        # Also dump the electrostatics
        if self.do_electrostatics:
            electrostatics_filename = os.path.join(self.path, self.electrostatics_name)
            self.gromacs_electrostatics.electrostatic_df.to_hdf(electrostatics_filename, key='electrostatics_df', mode='w')
            if self.verbose: print(self.gromacs_electrostatics.electrostatic_df)

        # Also also dump the forces
        if self.do_electrostatics:
            forces_filename = os.path.join(self.path, self.forces_name)
            self.master_forces_df.to_hdf(forces_filename, key='master_forces_df', mode='w')
            if self.verbose: print(self.master_forces_df)

        # Dump the pickle file(s)
        pickle_filename = os.path.join(self.path, self.pickle_name)
        with open(pickle_filename, 'wb') as f:
            pickle.dump(self.helix_analysis,f)

        if self.verbose: print(f"GromacsSeed::WriteData return")

    def WritePLUMEDAlphaRMSD(self):
        r""" Write the alpha_rmsd file for plumed to use
        """
        if self.verbose: print(f"GromacsSeed::WritePLUMEDAlphaRMSD")
        # Create the plumed file
        with open(os.path.join(self.path, self.plumed_file), "w") as fhandle:
            fstring = r'''# vim:ft=plumed
MOLINFO STRUCTURE={}

# Get detailed timing
DEBUG DETAILED_TIMERS

# Figure out the lipid and helix COM coordinates
# lipid: center of mass of phosphorus only
lipid_com: COM ATOMS={}
helix_com: COM ATOMS={}

# If we have a wholehelix, do that as well
wholehelix_com: COM ATOMS={}

# Get the Z distance
z_dist: DISTANCE ATOMS=lipid_com,helix_com COMPONENTS NOPBC
dz: COMBINE ARG=z_dist.z PERIODIC=NO

# Get the Z distance to the whole helix
z_dist_whole: DISTANCE ATOMS=lipid_com,wholehelix_com COMPONENTS NOPBC
wholedz: COMBINE ARG=z_dist_whole.z PERIODIC=NO

# Get the alpha value
alpha: ALPHARMSD RESIDUES={}

# Get the alpha values for the whole helix
wholealpha: ALPHARMSD RESIDUES={}

# Print to a file
PRINT ARG=dz,alpha,wholedz,wholealpha FILE={} STRIDE=1
'''.format(self.plumed_ref,
           ','.join(str(x) for x in self.lipid_P_atoms),
           ','.join(str(x) for x in self.helix_C_atoms),
           ','.join(str(x) for x in self.wholehelix_C_atoms),
           self.alpha_residues,
           self.wholehelix_residues,
           self.colvar_file)
            fhandle.write(fstring)

        if self.verbose: print(f"GromacsSeed::WritePLUMEDAlphaRMSD return")

    def GraphZdist(self, ax, color = 'b'):
        r""" Plot the absolute Z position of the helix
        """
        graph_seed_zdist(self, ax, color = color)

    def GraphZpos(self, ax, color = 'b'):
        r""" Plot the Z position of the lipid headgroups and protein
        """
        graph_seed_zpos_wheads(self, ax, color = color, docolvar = True)

    def GraphZcorr(self, ax, color = 'b'):
        r""" Plot the Z self-correlation (time lag) for a given deltaT
        """
        graph_seed_zcorr(self, ax, color = color)

    def GraphHelicity(self, ax, color = 'b'):
        r""" Plot the fractional helicity of the helix
        """
        graph_seed_helicity(self, ax, color = 'b', docolvar = True) 

    def GraphAlpha(self, ax, color = 'b'):
        r""" Plot the raw alpha (DSSP and RMSD) values for the helix
        """
        graph_seed_alpha(self, ax, color = 'b', docolvar = True)

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

    def GraphZForce(self, ax, color = 'b'):
        r""" Plot the Z forces on the helix (varying shades of alpha)
        """
        graph_seed_zforce(self, ax, color = 'b')

    def GraphPerpTorque(self, ax, color = 'b'):
        r""" Plot the perpendicular torque (from helix uhat and z-axis)
        """
        graph_seed_perptorque(self, ax, color = 'b')

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

