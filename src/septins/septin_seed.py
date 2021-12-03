# XXX: Put a license here

"""Class for a single septin seed"""

import os
import random
import sys
import yaml

import gsd.hoomd
import numpy as np
import pandas as pd

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
#from block_copolymer import BlockCopolymer
from membrane import Membrane
from block_ahdomain import BlockAHDomain
from helix_ahdomain import HelixAHDomain
from seed_base import SeedBase
from seed_graph_funcs import *

class SeptinSeed(SeedBase):
    def __init__(self, path, opts):
        SeedBase.__init__(self, path, opts)

        # Set the internal data classes
        self.lipids = None
        self.ahdomain = None

        self.ReadData()

        # Set the RNG state
        random.seed(self.nseed)

        # Lipid parameters
        if 'membrane' in self.default_yaml:
            self.lipids = Membrane(self.bead_size, self.default_yaml)
        # AH parameters
        if 'ah_domain' in self.default_yaml:
            self.ahtype = self.default_yaml['ah_domain']['polymer_type']
            if self.ahtype == 'block_copolymer':
                self.ahdomain = BlockAHDomain(self.bead_size, self.default_yaml)
            elif self.ahtype == 'helix_block_copolymer':
                self.ahdomain = HelixAHDomain(self.bead_size, self.default_yaml)
            else:
                print(f"AH domain type {self.ahtype} not implemented, exiting!")
                sys.exit(1)

        ## Set the trajectory file we are analyzing
        #self.trajectory_filename = self.default_yaml['simulation']['trajectory_file']

        ## Get membrane information
        #self.nbeads_membrane = np.float64(self.default_yaml['membrane']['nbeads'])
        #self.nlipids_per_layer = np.square(np.int64(self.default_yaml['membrane']['ngrid']))
        #self.nlipids = 2*self.nlipids_per_layer
        #self.nmembrane = self.nlipids * self.nbeads_membrane

        ## Get AH information if it exists
        #if 'ah_domain' in self.default_yaml:
        #    self.nahdomains = np.int64(self.default_yaml['ah_domain']['n_ah'])
        #    self.ah_start_nx = np.int64(self.nmembrane)
        #    self.ahdomain = BlockCopolymer(self.default_yaml, self.ah_start_nx)


        #self.hd5_name = "SeptinSeed.h5"

        #print(f"  Membrane beads: {self.nmembrane}")
        #print(f"  AH starting index: {self.ah_start_nx}")

    def ReadData(self):
        r""" Read the data from a YAML file for a single seed
        """
        # Simulation parameters
        self.kT                 = np.float64(self.default_yaml['simulation']['kT'])
        self.bead_size          = np.float64(self.default_yaml['simulation']['bead_size'])
        self.deltatau           = np.float64(self.default_yaml['simulation']['deltatau'])
        self.nsteps             = np.int32(np.float64(self.default_yaml['simulation']['nsteps']))
        self.nwrite             = np.int32(np.float64(self.default_yaml['simulation']['nwrite']))
        self.lbox               = np.float64(self.default_yaml['simulation']['lbox'])
        self.nseed              = np.int32(np.float64(self.default_yaml['simulation']['seed']))
        self.compute_mode       = self.default_yaml['simulation']['mode']
        self.trajectory_file    = self.default_yaml['simulation']['trajectory_file']
        self.init_type          = self.default_yaml['simulation']['init_type']
        self.integrator         = self.default_yaml['simulation']['integrator']

        # Check integrator settings to make sure it's okay!
        if (self.integrator != 'langevin' and self.integrator != 'NPT' and self.integrator != 'NPH'):
            print(f"Integrator {self.integrator} not used for this system, exiting!")
            sys.exit(1)

        # Unfortunately, now have some gymnastics for checking if something exists and setting
        if 'pdamp' in self.default_yaml['simulation']:
            self.pdamp = np.float64(self.default_yaml['simulation']['pdamp'])
        else:
            self.pdamp = None

        if 'tau' in self.default_yaml['simulation']:
            self.tau = np.float64(self.default_yaml['simulation']['tau'])
        else:
            self.tau = None

        if 'tauS' in self.default_yaml['simulation']:
            self.tauS = np.float64(self.default_yaml['simulation']['tauS'])
        else:
            self.tauS = None

        # Are we reading in previous information?
        if self.init_type == 'all':
            self.init_filename = ''
        elif self.init_type == 'read_gsd':
            self.init_filename = self.default_yaml['simulation']['init_filename']
        else:
            print(f"Need to specify a correct initialization type, tried {self.init_type}, exiting!")
            sys.exit(1)

    def PrintInformation(self, snap):
        r""" Print parameter information for SeptinSeed
        """
        print(f"--------")
        print(f"System information")
        print(f"Compute mode            = {self.compute_mode}")
        print(f"Integrator              = {self.integrator}")
        print(f"Delta tau               = {self.deltatau}")
        print(f"Simulation time (tau)   = {self.deltatau*self.nsteps}")
        print(f"kBT                     = {self.kT}")
        print(f"seed                    = {self.nseed}")
        print(f"box                     = ({snap.configuration.box[0]}, {snap.configuration.box[1]}, {snap.configuration.box[2]})")
        # Configurable parameters
        if self.tau: print(f"tau                     = {self.tau}")
        if self.pdamp: print(f"pdamp                   = {self.pdamp}")
        if self.tauS: print(f"tauS                    = {self.tauS}")

        if self.lipids:
            self.lipids.PrintInformation(snap)
        if self.ahdomain:
            self.ahdomain.PrintInformation(snap)

    def Configure(self, snap):
        r""" Configure the simulation snapshot for HOOMD
        """
        # Create a default configuration box size if we don't have a membrane
        snap.configuration.box = [self.lbox, self.lbox, self.lbox, 0, 0, 0]

        if self.lipids:
            self.lipids.InitMembrane(snap)
        if self.ahdomain:
            self.ahdomain.InitAH(snap)

    def GetHeadIndices(self, traj):
        r"""Get the head indices for the lipids, both upper
        and lower leaflets
        """
        head_idx = np.arange(0, self.nmembrane, self.nbeads)
        idx_skip = self.nbeads*2
        leaf1_idx = np.arange(0, self.nmembrane, idx_skip)
        lead2_idx = np.arange(self.nbeads, self.nmembrane, idx_skip)

        return [head_idx, leaf1_idx, leaf2_idx]

    def CheckLoadAnalyze(self, file_path, force_analyze = False):
        r""" Check if the data can be loaded from an HD5 file
        """
        print(f"Forcing analysis: {force_analyze}")
        if os.path.isfile(file_path) and not force_analyze:
            # Skip analysis and load
            try:
                self.LoadHD5(file_path)
                return False
            except EOFError: return False
            except: raise
        else:
            return True

    def SaveHD5(self):
        r""" Save the current state in HD5 format
        """
        hd5_filename = os.path.join(self.path, "SeptinSeed.h5")
        self.df.to_hdf(hd5_filename, key='df', mode='w')

    def LoadHD5(self, file_path):
        r""" Load the current state in HD5 format
        """
        self.df = pd.read_hdf(file_path)

    def Analyze(self, force_analyze = False):
        r"""Analysis of septin/membrane simulation seed
        """
        # Check the options for if we are forcing analysis
        if self.opts.force:
            force_analyze = True
        # Check if we are loading or analyzing the information
        if not self.CheckLoadAnalyze(os.path.join(self.path, self.hd5_name),
                                     force_analyze):
            print(f"Loading seed information")
            return False

        print(f"Running seed analysis")
        print(f"Path: {self.path}")
        print(f"  yaml file: {self.yaml_filename}")
        print(f"  trajectory file: {self.trajectory_filename}")

        # If we are doing the analysis load the trajectory file
        self.traj_all = gsd.hoomd.open(os.path.join(self.path, self.trajectory_filename))
        print(f"Analyzing nframes: {len(self.traj_all)}")

        # XXX: These are hardcoded for now, change later for the analysis for what frames, etc
        # we want to do the analysis on
        self.nmeasure = 10
        self.min_frame = 000
        self.max_frame = 200

        # Set up the storage arrays for variables
        self.timedata = {}
        self.timedata['timestep'] = [] # Create an empty list for the frame times


        # Main analysis loop
        for itx,traj in enumerate(self.traj_all):
            # Bail if we don't want to analyze this timepoint
            if itx % self.nmeasure != 0:
                continue
            if itx < self.min_frame:
                continue
            if itx >= self.max_frame:
                break
            print(f"    Analyzing frame: {itx}")

            # Get the timestep of the current analysis from the trajectory
            timestep = np.int64(traj.configuration.step)
            print(f"    Timestep: {timestep}")
            self.timedata['timestep'].append(timestep)

            # Here are the actual guts of the analyses
            self.Temperature(timestep, traj)
            self.Pressure(timestep, traj)
            self.MembraneArea(timestep, traj)
            




        # Make sure to clean up after ourselves if we've opened the file
        self.traj_all.close()

        # Now do the postprocessing and convert everything into a pandas dataframe
        self.Postprocess()

        # Save the data
        self.SaveHD5()

    def Postprocess(self):
        r""" Postprocess the results to save in HD5F format
        """
        # Have to do timestep separately
        timestep = np.array(self.timedata['timestep'], dtype = np.int64)
        df_timestep = pd.DataFrame(timestep, columns = ['timestep'])

        dfs = []
        dfs.append(df_timestep)
        for key in sorted(self.timedata):
            if key == 'timestep':
                continue
            val_dict = self.timedata[key]
            val_arr = np.array([val_dict[ts] for ts in timestep], dtype = np.float64)
            df = pd.DataFrame(val_arr, columns = [key])
            dfs.append(df)

        self.df = pd.concat(dfs, axis=1)
        print(self.df)

    def Temperature(self, timestep, traj):
        r""" Simulation kinetic temperature
        """
        T = traj.log['md/compute/ThermodynamicQuantities/kinetic_temperature'][0]
        if 'T' not in self.timedata:
            self.timedata['T'] = {}
        self.timedata['T'][timestep] = T

    def Pressure(self, timestep, traj):
        r""" Simulation pressure
        """
        P = traj.log['md/compute/ThermodynamicQuantities/pressure'][0]
        if 'P' not in self.timedata:
            self.timedata['P'] = {}
        self.timedata['P'][timestep] = P

    def MembraneArea(self, timestep, traj):
        r""" Membrane area
        """
        Lx = np.float64(traj.configuration.box[0])
        Ly = np.float64(traj.configuration.box[1])
        membrane_area = Lx*Ly
        if 'membrane_area' not in self.timedata:
            self.timedata['membrane_area'] = {}
        self.timedata['membrane_area'][timestep] = membrane_area

    def Graph(self, axarr, color = 'b'):
        r"""Default graphing call for single seeds
        """
        graph_scatter(self.label, self.df["timestep"], self.df["T"], axarr[0], mtitle = "Temperature", xtitle = "Timestep", ytitle = "Temperature (kT)")
        graph_scatter(self.label, self.df["timestep"], self.df["P"], axarr[1], mtitle = "Pressure", xtitle = "Timestep", ytitle = "Presure (kT/$\sigma^{3}$)")
        graph_scatter(self.label, self.df["timestep"], self.df["membrane_area"], axarr[2], mtitle = "Membrane Area", xtitle = "Timestep", ytitle = "Membrane area ($\sigma^{2}$)")


