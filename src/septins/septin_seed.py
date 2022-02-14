# XXX: Put a license here

"""Class for a single septin seed"""

import os
import random
import sys
import yaml

import gsd.hoomd
import numpy as np
import pandas as pd

# Analysis imports
import freud
from freud import box

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from simple_spheres import SimpleSpheres
from membrane import Membrane
from block_ahdomain import BindingState, BlockAHDomain
from helix_ahdomain import HelixAHDomain
from seed_base import SeedBase
from common import radial_average
from seed_graph_funcs import *

class SeptinSeed(SeedBase):
    def __init__(self, path, opts):
        print(f"SeptinSeed::__init__")

        SeedBase.__init__(self, path, opts)

        # Set the internal data classes
        self.lipids = None
        self.ahdomain = None
        self.simple_spheres = None

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
        # Simple spheres for testing purposes
        if 'simple_spheres' in self.default_yaml:
            self.simple_spheres = SimpleSpheres(self.bead_size, self.default_yaml)

        # Set the name for the results file at the end
        self.hd5_name = "SeptinSeed.h5"

        print(f"SeptinSeed::__init__ return")

    def ReadData(self):
        r""" Read the data from a YAML file for a single seed
        """
        print(f"SeptinSeed::ReadData")
        # Simulation parameters
        self.kT                 = np.float64(self.default_yaml['simulation']['kT'])
        self.bead_size          = np.float64(self.default_yaml['simulation']['bead_size'])
        self.deltatau           = np.float64(self.default_yaml['simulation']['deltatau'])
        self.nsteps             = np.int32(np.float64(self.default_yaml['simulation']['nsteps']))
        self.nwrite             = np.int32(np.float64(self.default_yaml['simulation']['nwrite']))
        self.min_frame          = np.int32(np.float64(self.default_yaml['simulation']['min_frame']))
        self.max_frame          = np.int32(np.float64(self.default_yaml['simulation']['max_frame']))
        self.lbox               = np.float64(self.default_yaml['simulation']['lbox'])
        self.nseed              = np.int32(np.float64(self.default_yaml['simulation']['seed']))
        self.compute_mode       = self.default_yaml['simulation']['mode']
        self.trajectory_file    = self.default_yaml['simulation']['trajectory_file']
        self.init_type          = self.default_yaml['simulation']['init_type']
        self.integrator         = self.default_yaml['simulation']['integrator']
        self.time_divisions     = np.int32(np.float64(self.default_yaml['simulation']['time_divisions']))

        # Check integrator settings to make sure it's okay!
        if (self.integrator != 'langevin' and self.integrator != 'NPT' and self.integrator != 'NPH' and self.integrator != 'Brownian'):
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
        print(f"SeptinSeed::ReadData return")

    def PrintInformation(self, snap):
        r""" Print parameter information for SeptinSeed
        """
        print(f"--------")
        print(f"System information")
        print(f"Compute mode            = {self.compute_mode}")
        print(f"Integrator              = {self.integrator}")
        print(f"Delta tau               = {self.deltatau}")
        print(f"Simulation time (tau)   = {self.deltatau*self.nsteps}")
        print(f"nwrite                  = {self.nwrite}")
        print(f"min/max frame           = {self.min_frame}, {self.max_frame}")
        print(f"time divisions          = {self.time_divisions}")
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
        if self.simple_spheres:
            self.simple_spheres.PrintInformation(snap)

    def Configure(self, snap):
        r""" Configure the simulation snapshot for HOOMD
        """
        print(f"SeptinSeed::Configure")
        # Create a default configuration box size if we don't have a membrane
        snap.configuration.box = [self.lbox, self.lbox, self.lbox, 0, 0, 0]

        if self.lipids:
            self.lipids.InitMembrane(snap)
        if self.ahdomain:
            self.ahdomain.InitAH(snap)
        if self.simple_spheres:
            self.simple_spheres.InitSimpleSpheres(snap)

        print(f"SeptinSeed::Configure return")

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
        print(f"SeptinSeed::Analyze")
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
        print(f"  trajectory file: {self.trajectory_file}")

        # If we are doing the analysis load the trajectory file
        self.traj_all = gsd.hoomd.open(os.path.join(self.path, self.trajectory_file))
        print(f"Analyzing nframes: {len(self.traj_all)}")

        # Configure membrane/ahdomain and print the information
        if self.lipids:
            self.lipids.ConfigureAnalysis(self.traj_all[0])
        if self.ahdomain:
            self.ahdomain.ConfigureAnalysis(self.traj_all[0])
        if self.simple_spheres:
            self.simple_spheres.ConfigureAnalysis(self.traj_all[0])
        self.PrintInformation(self.traj_all[0])

        # Set up the storage arrays for variables
        self.timedata = {}
        self.modedata = {}
        self.distdata = {}
        self.timedata['timestep'] = [] # Create an empty list for the frame times

        # Keep track of the current frame and the current time division
        self.current_time_division = 0
        self.current_time_trigger = False
        self.first_trigger = True

        # Main analysis loop
        for itx,snap in enumerate(self.traj_all):
            # Analyze every timepoint in the sequence
            self.current_frame =  itx # Keep track of the current frame number
            #print(f"division: {(self.max_frame - self.min_frame)/self.time_divisions}")
            if (self.current_frame % np.int64((self.max_frame - self.min_frame)/self.time_divisions) == 0 and not self.first_trigger):
                self.current_time_trigger = True # XXX Probably a better way to do this, but need to get it done
            elif (self.first_trigger):
                self.first_trigger = False

            if itx < self.min_frame:
                continue
            if itx >= self.max_frame:
                break
            print(f"Analyzing frame: {itx}")

            # Get the timestep of the current analysis from the trajectory
            timestep = np.int64(snap.configuration.step)
            if self.current_time_trigger:
                print(f"    Timestep: {timestep}, time_trigger: {self.current_time_trigger}, {self.current_time_division}")
            else:
                print(f"    Timestep: {timestep}")

            self.timedata['timestep'].append(timestep)

            # Here are the actual guts of the analyses
            self.Temperature(timestep, snap)
            self.Pressure(timestep, snap)
            self.MembraneArea(timestep, snap)

            # Membrane analysis
            self.MembraneModes(timestep, snap)

            # AH domain analysis
            self.AHDomainKon(timestep, snap)

            # Simple spheres analysis
            self.SimpleDiffusion(timestep, snap)

            # Reset the triggers
            if self.current_time_trigger:
                self.current_time_division += 1
            self.current_time_trigger = False

        # Make sure to clean up after ourselves if we've opened the file
        self.traj_all.close()

        # Now do the postprocessing and convert everything into a pandas dataframe
        self.Postprocess()

        # Save the data
        self.SaveHD5()

        print(f"SeptinSeed::Analyze return")

    def Postprocess(self):
        r""" Postprocess information for things like distributions
        """

        # Increment the time_division counter to force a change where needed
        self.current_time_division += 1

        # Postprocess dynamic data (e.g. MSD for simple spheres)
        self.PostprocessDynamics()

        # Postprocess distributions
        self.PostprocessDistributions()

        # Also postprocess for saving
        self.PostprocessHD5()

    def PostprocessDynamics(self):
        r""" Postprocess dynamic data
        """
        # Early return if something doesn't exist
        if 'simple_diffusion' not in self.timedata:
            return

        timestep = np.array(self.timedata['timestep'], dtype = np.int64)
        positions = self.timedata['simple_diffusion']
        positions_arr = np.array([positions[ts] for ts in timestep], dtype = np.float64)

        print(f"{positions_arr.shape}")

        final_box = freud.Box.cube(100.0)
        msd_compute = freud.msd.MSD(box = final_box, mode = 'window')
        msd_compute.compute(positions_arr)
        msd_arr = msd_compute.msd

        # Enter this into the timedata calculation
        self.timedata['msd_simple'] = {}
        for idx,ts in enumerate(timestep):
            self.timedata['msd_simple'][ts] = msd_arr[idx]

    def PostprocessHD5(self):
        r""" Postprocess the results to save in HD5F format
        """
        # Have to do timestep separately
        timestep = np.array(self.timedata['timestep'], dtype = np.int64)
        df_timestep = pd.DataFrame(timestep, columns = ['timestep'])

        # Process the timesteps
        dfs = []
        dfs.append(df_timestep)
        for key in sorted(self.timedata):
            if key == 'timestep':
                continue
            if key == 'simple_diffusion':
                continue
            val_dict = self.timedata[key]
            val_arr = np.array([val_dict[ts] for ts in timestep], dtype = np.float64)
            df = pd.DataFrame(val_arr, columns = [key])
            dfs.append(df)

        # Process the distribution data too
        if self.lipids and self.ahdomain:
            dfs.append(self.df_lifetimes)
            dfs.append(self.subdivision_df)

        # Process membrane mode data
        if self.lipids:
            dfs.append(self.df_x_fft)
            dfs.append(self.df_su_fft)
            dfs.append(self.df_qcutoff_fft)
            dfs.append(self.df_x_direct)
            dfs.append(self.df_su_direct)
            dfs.append(self.df_qcutoff_direct)

        self.df = pd.concat(dfs, axis=1)
        print(self.df)

    def PostprocessDistributions(self):
        r""" Postprocess the distributions
        """
        self.PostprocessLifetimes()
        self.PostprocessMembraneModes()

    def PostprocessLifetimes(self):
        r""" Postprocess lifetime distributions
        """
        # Early exit if it doesn't exist
        if not self.lipids or not self.ahdomain:
            return

        # Postprocess lifetimes
        lifetime_dict = {}
        lifetime_dict[BindingState.free] = []
        lifetime_dict[BindingState.near] = []
        lifetime_dict[BindingState.surface] = []
        lifetime_dict[BindingState.intermediate] = []
        lifetime_dict[BindingState.deep] = []
        subdivision_lifetime_dict = {}
        # Do the lifetime distribution
        for ahdx in range(self.ahdomain.nah):
        #for ahdx in range(1):
            subdomain = self.ahdomain.GetSubAHDomain(ahdx)
            # Flush the binding state information
            subdomain.UpdateBindingState(-1, self.current_time_division)

            lifetime_dict[BindingState.free].extend(subdomain.lifetime_dict[BindingState.free]) 
            lifetime_dict[BindingState.near].extend(subdomain.lifetime_dict[BindingState.near]) 
            lifetime_dict[BindingState.surface].extend(subdomain.lifetime_dict[BindingState.surface]) 
            lifetime_dict[BindingState.intermediate].extend(subdomain.lifetime_dict[BindingState.intermediate]) 
            lifetime_dict[BindingState.deep].extend(subdomain.lifetime_dict[BindingState.deep]) 

            # Now compile the results for all the subdivision, the first division counts at the end of the range
            for itx in range(0, self.current_time_division):
                if itx not in subdivision_lifetime_dict:
                    subdivision_lifetime_dict[itx] = {}
                    subdivision_lifetime_dict[itx][BindingState.free] = []
                    subdivision_lifetime_dict[itx][BindingState.near] = []
                    subdivision_lifetime_dict[itx][BindingState.surface] = []
                    subdivision_lifetime_dict[itx][BindingState.intermediate] = []
                    subdivision_lifetime_dict[itx][BindingState.deep] = []
                subdivision_lifetime_dict[itx][BindingState.free].extend(subdomain.lifetime_subdivision_dict[itx][BindingState.free]) 
                subdivision_lifetime_dict[itx][BindingState.near].extend(subdomain.lifetime_subdivision_dict[itx][BindingState.near]) 
                subdivision_lifetime_dict[itx][BindingState.surface].extend(subdomain.lifetime_subdivision_dict[itx][BindingState.surface]) 
                subdivision_lifetime_dict[itx][BindingState.intermediate].extend(subdomain.lifetime_subdivision_dict[itx][BindingState.intermediate]) 
                subdivision_lifetime_dict[itx][BindingState.deep].extend(subdomain.lifetime_subdivision_dict[itx][BindingState.deep]) 

        # Convert lifetime information to DFs to save
        self.df_lifetimes = pd.DataFrame(dict([ (k.name,pd.Series(np.float32(v))) for k,v in lifetime_dict.items() ]))

        # Conver the time subdivisions into DFs to save
        subdivision_dfs = []
        for time_division,bstate in subdivision_lifetime_dict.items():
            # For each state name, create a new DF
            for state_name,val in bstate.items():
                col_name = state_name.name + str(time_division)
                df = pd.DataFrame(val, columns = [col_name])
                subdivision_dfs.append(df)
        self.subdivision_df = pd.concat(subdivision_dfs, axis=1)

    def PostprocessMembraneModes(self):
        r""" Postprocess membrane mode data
        """
        # Early check if we have somehow not run the analysis
        if not self.lipids:
            return

        # Set up the bin size
        deltaq = 0.0375 # Should be 0.05 nm^-1

        # Loop over the timesteps
        timestep = np.array(self.timedata['timestep'], dtype=np.int64)

        membranemodes_fft = self.modedata['membranemodes_fft']
        membranemodes_fft_arr = np.array([membranemodes_fft[ts] for ts in timestep], dtype = np.complex128)
        qcutoff_fft = self.modedata['qcutoff_fft']
        qcutoff_fft_arr = np.array([qcutoff_fft[ts] for ts in timestep], dtype = np.float64)
        membranemodes_direct = self.modedata['membranemodes_direct']
        membranemodes_direct_arr = np.array([membranemodes_direct[ts] for ts in timestep], dtype = np.complex128)
        qcutoff_direct = self.modedata['qcutoff_direct']
        qcutoff_direct_arr = np.array([qcutoff_direct[ts] for ts in timestep], dtype = np.float64)

        # Loop over the membrane modes
        intensity_fft_list = []
        intensity_direct_list = []
        qcutoff_fft_list = []
        qcutoff_direct_list = []
        for itx in np.arange(membranemodes_fft_arr.shape[0]):
            [radii_fft, intensity_fft] = radial_average(membranemodes_fft_arr[itx,:,:], deltaq, qcutoff_fft_arr[itx])
            intensity_fft_list.append(intensity_fft)
            qcutoff_fft_list.append(qcutoff_fft_arr[itx])
            
            [radii_direct, intensity_direct] = radial_average(membranemodes_direct_arr[itx,:,:], deltaq, qcutoff_direct_arr[itx])
            intensity_direct_list.append(intensity_direct)
            qcutoff_direct_list.append(qcutoff_direct_arr[itx])

        uq_fft_arr      = np.array(intensity_fft_list)
        uq_direct_arr   = np.array(intensity_direct_list)

        uq_fft_mean     = np.mean(uq_fft_arr, axis = 0)
        uq_direct_mean  = np.mean(uq_direct_arr, axis = 0)

        su_fft          = np.square(uq_fft_mean)*self.lipids.nlipids_per_leaflet
        su_direct       = np.square(uq_direct_mean)*self.lipids.nlipids_per_leaflet

        self.df_x_fft           = pd.DataFrame(radii_fft, columns=['x_fft'])
        self.df_su_fft          = pd.DataFrame(su_fft, columns=['su_fft'])
        self.df_qcutoff_fft     = pd.DataFrame(qcutoff_fft_list, columns=['qcutoff_fft'])
        self.df_x_direct        = pd.DataFrame(radii_direct, columns=['x_direct'])
        self.df_su_direct       = pd.DataFrame(su_direct, columns=['su_direct'])
        self.df_qcutoff_direct  = pd.DataFrame(qcutoff_direct_list, columns=['qcutoff_direct'])

    def Temperature(self, timestep, snap):
        r""" Simulation kinetic temperature
        """
        T = snap.log['md/compute/ThermodynamicQuantities/kinetic_temperature'][0]
        if 'T' not in self.timedata:
            self.timedata['T'] = {}
        self.timedata['T'][timestep] = T

    def Pressure(self, timestep, snap):
        r""" Simulation pressure
        """
        P = snap.log['md/compute/ThermodynamicQuantities/pressure'][0]
        if 'P' not in self.timedata:
            self.timedata['P'] = {}
        self.timedata['P'][timestep] = P

    def MembraneArea(self, timestep, snap):
        r""" Membrane area
        """
        Lx = np.float64(snap.configuration.box[0])
        Ly = np.float64(snap.configuration.box[1])
        membrane_area = Lx*Ly
        if 'membrane_area' not in self.timedata:
            self.timedata['membrane_area'] = {}
        self.timedata['membrane_area'][timestep] = membrane_area

    def SimpleDiffusion(self, timestep, snap):
        r""" Diffusion analysis of parts of the system
        """

        if not self.simple_spheres:
            return

        if 'simple_diffusion' not in self.timedata:
            self.timedata['simple_diffusion'] = {}
        s_idx = self.simple_spheres.GetSimpleSpheresIndices(snap)
        self.timedata['simple_diffusion'][timestep] = snap.particles.position[s_idx]

    def MembraneModes(self, timestep, snap):
        r""" Membrane bending modes for both direct and fft measurements
        """

        if not self.lipids:
            return

        # Constants used for lots of things
        Nx = 100
        Ny = 100
        Nxgrid = Nx * 1j
        Nygrid = Ny * 1j
        Ndirect = 5 # Number of direct measurements to make

        # Get the head, intermeidate, and tail indices for the lipids
        [h_idx, leaf1_h_idx, leaf2_h_idx] = self.lipids.GetLeafletIndices(snap, 'H')
        # Get the box size for this timeframe
        Lx = np.float64(snap.configuration.box[0])
        Ly = np.float64(snap.configuration.box[1])
        qcutoffx = 2.0*np.pi/Lx
        qcutoffy = 2.0*np.pi/Ly

        ################
        # Interpolation method on a grid
        ################
        # Get the Z positions of the lipid heads and recenter
        positions = snap.particles.position
        z1 = positions[leaf1_h_idx, 2]
        z2 = positions[leaf2_h_idx, 2]
        z0 = (np.sum(z1) + np.sum(z2))/(len(z1) + len(z2))
        z1 = z1 - z0
        z2 = z2 - z0

        # Get the r positions for each leaflet
        r1 = positions[leaf1_h_idx, 0:2]
        r2 = positions[leaf2_h_idx, 0:2]

        # Interpolation method to put on a grid
        # XXX: Hardcoded for now, change later!
        xpixel = Lx/Nx
        ypixel = Ly/Ny
        grid_x, grid_y = np.mgrid[-Lx/2:Lx/2:Nxgrid, -Ly/2:Ly/2:Nygrid]
        from scipy.interpolate import griddata
        grid_z1 = griddata(r1, z1, (grid_x, grid_y), method = 'cubic', fill_value = np.mean(z1))
        grid_z2 = griddata(r2, z2, (grid_x, grid_y), method = 'cubic', fill_value = np.mean(z2))

        # Take the 2d fourier transform
        u1 = np.fft.fft2(grid_z1, norm = 'forward')
        u1shift = np.fft.fftshift(u1)
        u2 = np.fft.fft2(grid_z2, norm = 'forward')
        u2shift = np.fft.fftshift(u2)

        # Get the sample frequencies
        freqx = np.fft.fftshift(np.fft.fftfreq(u1shift.shape[1], xpixel))
        freqy = np.fft.fftshift(np.fft.fftfreq(u2shift.shape[0], ypixel))

        # Compute final variables
        uq_2d_fft = 0.5*(u1shift + u2shift)

        ################
        # Direct measurement of fourier coefficients
        ################
        u1direct = np.zeros((2*Ndirect+1,2*Ndirect+1), dtype=np.complex128)
        u2direct = np.zeros((2*Ndirect+1,2*Ndirect+1), dtype=np.complex128)

        for n in np.arange(-Ndirect,Ndirect+1,1):
            for m in np.arange(-Ndirect,Ndirect+1,1):
                idx = n + Ndirect
                kdx = m + Ndirect

                q = 2.0*np.pi*np.array([n/Lx, m/Ly])

                for k,r1k in enumerate(r1):
                    val = z1[k] * np.exp(-1j*np.dot(q,r1k))
                    u1direct[idx,kdx] += val
                for k,r2k in enumerate(r2):
                    val = z2[k] * np.exp(-1j*np.dot(q,r2k))
                    u2direct[idx,kdx] += val

        # Compute the final version
        uq_2d_direct = 1.0/(2.0*self.lipids.nlipids_per_leaflet)*(u1direct + u2direct)

        # Save the information we need later
        if 'membranemodes_fft' not in self.modedata:
            self.modedata['membranemodes_fft'] = {}
        if 'membranemodes_direct' not in self.modedata:
            self.modedata['membranemodes_direct'] = {}
        if 'qcutoff_fft' not in self.modedata:
            self.modedata['qcutoff_fft'] = {}
        if 'qcutoff_direct' not in self.modedata:
            self.modedata['qcutoff_direct'] = {}

        # For the FFT, make sure to convert to q-space (has the 2pi in the numerator)
        self.modedata['membranemodes_fft'][timestep] = uq_2d_fft
        self.modedata['qcutoff_fft'][timestep] = (freqx[1] - freqx[0])*2.0*np.pi
        
        self.modedata['membranemodes_direct'][timestep] = uq_2d_direct
        self.modedata['qcutoff_direct'][timestep] = qcutoffx

    def AHDomainKon(self, timestep, snap):
        r""" k_on measurements for ah domains
        """
        # Bail if there isn't an AH domain or a membrane
        if not self.ahdomain or not self.lipids:
            return

        # Create the lifetime subdivisions if we need to
        if 'lifetime_subdivisions' not in self.distdata:
            self.distdata['lifetime_subdivisions'] = {}

        # Get the head, intermeidate, and tail indices for the lipids
        [h_idx, leaf1_h_idx, leaf2_h_idx] = self.lipids.GetLeafletIndices(snap, 'H')
        [i_idx, leaf1_i_idx, leaf2_i_idx] = self.lipids.GetLeafletIndices(snap, 'I')
        [t_idx, leaf1_t_idx, leaf2_t_idx] = self.lipids.GetLeafletIndices(snap, 'T')

        # Loop over membrane/AH combinations to look for the interactions,
        # and then record them in the individual SingleAHDomains
        for ahdx in range(self.ahdomain.nah):
        #for ahdx in range(1):
            # Get the current box size
            box_data = snap.configuration.box
            f_box = box.Box(Lx = box_data[0], Ly = box_data[1], Lz = box_data[2])

            subdomain = self.ahdomain.GetSubAHDomain(ahdx)

            # Set control variables first
            near_membrane = False
            surface_interaction = False
            intermediate_interaction = False
            deep_interaction = False

            r_cutoff = 1.0*(self.lipids.r0 + self.ahdomain.r0)

            # Check if we are near the membrane
            # XXX: There is probably a way to significantly speed this up by only constructing the
            #      AABB query once, and then filtering out the results?
            # XXX: For now, use a cutoff of something like within the range of 10 cutoff radii
            system_near = freud.AABBQuery(f_box, snap.particles.position[h_idx])
            query_points = snap.particles.position[subdomain.index_all]
            nlist_near = system_near.query(query_points, {"r_max": 10.0*r_cutoff}).toNeighborList()
            if len(nlist_near) > 0: near_membrane = True

            # Do the heads first, do we have interactions with the hydrophilic heads
            membrane_head_system = freud.AABBQuery(f_box, snap.particles.position[h_idx])
            # Check the hydrophilic portinos of the protein
            query_points = snap.particles.position[subdomain.index_ah1]
            nlist_surface = membrane_head_system.query(query_points, {"r_max": r_cutoff}).toNeighborList()
            if len(nlist_surface) > 0: surface_interaction = True

            # Check if hydrophobic regions are near tails
            membrane_intermediate_system = freud.AABBQuery(f_box, snap.particles.position[i_idx])
            query_points = snap.particles.position[subdomain.index_ah2]
            nlist_intermediate = membrane_intermediate_system.query(query_points, {"r_max": r_cutoff}).toNeighborList()
            if len(nlist_intermediate) > 0: intermediate_interaction = True

            # Check if hydrophobic regions are near the deep
            membrane_deep_system = freud.AABBQuery(f_box, snap.particles.position[t_idx])
            query_points = snap.particles.position[subdomain.index_ah2]
            nlist_deep = membrane_deep_system.query(query_points, {"r_max": r_cutoff}).toNeighborList()
            if len(nlist_deep) > 0: deep_interaction = True

            # Check all of the conditions, in order, to tell the ahdomain what it is doing
            # Set the binding type to what we see here in order
            # -1 for uninitialized, 0 for free,
            # 1 for near surface, 2 for surface interact,
            # 3 for intermediate, 4 for deep
            current_binding_state = BindingState.init
            if deep_interaction:
                current_binding_state = BindingState.deep
            elif intermediate_interaction:
                current_binding_state = BindingState.intermediate
            elif surface_interaction:
                current_binding_state = BindingState.surface
            elif near_membrane:
                current_binding_state = BindingState.near
            else:
                current_binding_state = BindingState.free

            # Now, pass this information to the subdomain to do with it what it will
            subdomain.UpdateBindingState(current_binding_state, self.current_time_division)

    def GraphDynamic(self, axarr, color = 'b'):
        r"""Default graphing call for single seeds dynamic information
        """
        graph_seed_scatter(self.label, self.df["timestep"], self.df["T"], axarr[0], mtitle = "Temperature", xtitle = "Timestep", ytitle = "Temperature (kT)")
        graph_seed_scatter(self.label, self.df["timestep"], self.df["P"], axarr[1], mtitle = "Pressure", xtitle = "Timestep", ytitle = "Presure (kT/$\sigma^{3}$)")
        graph_seed_scatter(self.label, self.df["timestep"], self.df["membrane_area"], axarr[2], mtitle = "Membrane Area", xtitle = "Timestep", ytitle = "Membrane area ($\sigma^{2}$)")

    def GraphDistributions(self, axarr, color = 'b'):
        r"""Default graphing call for single seed distribution information
        """
        if not self.lipids or not self.ahdomain:
            return
        graph_seed_lifetime_distribution(self.label, self.df[['free', 'near', 'surface', 'intermediate', 'deep']], axarr, mtitle = "Lifetime", xtitle = "Lifetime (frames)", ytitle = "Number")

    def GraphMembraneModes(self, ax, color = 'b'):
        r""" Default graphing acll for single seed membrane mode data
        """
        if not self.lipids:
            return
        graph_seed_membranemodes(self.label, self.df[['x_fft', 'su_fft', 'x_direct', 'su_direct', 'qcutoff_fft', 'membrane_area']], self.lipids.nlipids_per_leaflet, ax, mtitle = "Membrane Modes", xtitle = r'q ($\sigma^{-1}$)', ytitle = r'$ \langle | u(q) |^{2} \rangle $ ($\sigma^{2}$)')

    def GraphSimpleMSD(self, ax, color = 'b'):
        r""" Default graphing call for single particle MSD tests
        """
        if not self.simple_spheres:
            return
        graph_seed_simplespheres_msd(self.label, self.deltatau, self.df[['timestep', 'msd_simple']], ax, mtitle = "Simple MSD", xtitle = r'Timestep', ytitle = r'MSD ($\sigma^{2}$)')


