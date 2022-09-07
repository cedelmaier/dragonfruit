#!/usr/bin/env python3

# XXX: Put a license here

""" Simple script to combine the gromacs analyses listed, keep updating and clean up in future """

import pickle
import os
import sys

import MDAnalysis as mda
import MDAnalysis.transformations as trans
from MDAnalysis.analysis import helix_analysis as hel

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src', 'lib'))
from stylelib.common_styles import *

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data')
simulation_names = [
                    #'unbiased_zdepth00_rotx0_helix_50mMKCl',
                    #'unbiased_zdepth00_rotx90_helix_50mMKCl',
                    #'unbiased_zdepth00_rotx180_helix_50mMKCl',
                    #'unbiased_zdepth00_rotx270_helix_50mMKCl',
                    #'unbiased_zdepth00_trans_helix_50mMKCl',
                    #'unbiased_150mMKCl',
                    'rfmonomer_aglipid_11x11_zdepth00_rotx0_50mMKCl',
                    'rfmonomer_aglipid_11x11_zdepth00_rotx90_50mMKCl',
                    'rfmonomer_aglipid_11x11_zdepth00_rotx270_50mMKCl',
                    ]
simulation_legends = [
                      #r'Rot 0 deg',
                      #r'Rot 90 deg',
                      #r'Rot 180 deg',
                      #r'Rot 270 deg',
                      #r'Trans',
                      #r'Rot 0 deg 150 mM KCl',
                      r'RFMAGL 0 deg 50 KCl',
                      r'RFMAGL 90 deg 50 KCl',
                      r'RFMAGL 270 deg 50 KCl',
                     ]


# Set up an output name
#outname = 'coiled_all'
#outname = 'unfolded_main'
#outname = 'unfolded_10'
#outname = 'unfolded_20'
#outname = 'unbiased_zdepth00_rotxN_helix_50mMKCl' 
#outname = 'unbiased_zdepth00_rotxNtrans_helix_50mMKCl' 
#outname = 'unbiased_zdepth00_rotx90_helix_150mMKCl'
#outname = 'unbiased_zdepth00_rotx0_helix_50mMKCl'
#outname = 'unbiased_zdepth00_rotx90_helix_50mMKCl'
#outname = 'unbiased_zdepth00_rotx180_helix_50mMKCl'
#outname = 'unbiased_zdepth00_rotx270_helix_50mMKCl'
outname = 'rfmonomer_aglipid_rotxN_50mMKCl'

# Set up the plots beforehand
#plt.style.use(septin_poster_stl)
plt.style.use(septin_runs_stl)
fig_zdist, ax_zdist = plt.subplots(1, 1, figsize = (15, 10))
fig_zpos , ax_zpos  = plt.subplots(1, 1, figsize = (15, 10))
fig_twist, ax_twist = plt.subplots(1, 1, figsize = (15, 10))
fig_resid, ax_resid = plt.subplots(1, 1, figsize = (15, 10))
fig_helix, ax_helix = plt.subplots(1, 1, figsize = (15, 10))
fig_tilt , ax_tilt  = plt.subplots(1, 1, figsize = (15, 10))

# Set up a stride and only graph every stride points
stride = 10

# For each simulation, loop through and get the relevant data and plot it
Nsims = 4
cidx = 0
for isim,sname in enumerate(simulation_names):
    for iseed in np.arange(1, Nsims+1):
        print(f"Sim: {isim} = {sname}, N{iseed}")
        filepath = os.path.join(main_path, sname, f"N{iseed}")

        hd5_filename = os.path.join(filepath, f"N{iseed}.h5")
        master_time_df = pd.read_hdf(hd5_filename)

        times = master_time_df.index

        labelname = simulation_legends[isim] + f" N{iseed}"

        # Compute the zdist plot
        zdist = np.abs(master_time_df['helix_z'] - master_time_df['lipid_z'])
        ax_zdist.plot(times, zdist, label = labelname)

        # Compute the z position
        z_pos           = master_time_df[['helix_z']].to_numpy().flatten()
        leaflet0_pos    = master_time_df[['leaflet0_z']].to_numpy().flatten()
        leaflet1_pos    = master_time_df[['leaflet1_z']].to_numpy().flatten()
        lipid_pos       = master_time_df[['lipid_z']].to_numpy().flatten()
        z_pbc           = master_time_df[['unit_cell_z']].to_numpy().flatten()
        # Subtract off the lipid COM position
        z_pos = z_pos - lipid_pos
        leaflet0_pos = leaflet0_pos - lipid_pos
        leaflet1_pos = leaflet1_pos - lipid_pos
        z_pbc = z_pbc/2 - lipid_pos
        # Correct the position if under the lower leaflet
        for idx in range(len(z_pos)):
            if z_pos[idx] < leaflet1_pos[idx]:
                z_pos[idx] = z_pbc[idx] - z_pos[idx]

        # Plot the location of everything at some stride distance
        ax_zpos.plot(times[::stride], z_pos[::stride], label = labelname)
        ax_zpos.plot(times[::stride], leaflet0_pos[::stride], color = 'k')
        ax_zpos.plot(times[::stride], leaflet1_pos[::stride], color = 'k')

        # Compute the helicity
        helicity = master_time_df[['helicity']].to_numpy().flatten()
        # Smooth this data, as it is gross right now
        ax_helix.plot(times[::stride], helicity[::stride], label = labelname)

        # Get the tilt
        global_tilt = master_time_df[['global_tilt']].to_numpy().flatten()
        ax_tilt.plot(times[::stride], global_tilt[::stride], label = labelname)


    cidx += 1

ax_zdist.set_xlabel('Time (ps)')
ax_zdist.set_ylabel(r'Z distance ($\AA$)')
#ax_zdist.legend(loc = 'lower right')
fig_zdist.tight_layout()
fig_zdist.savefig('gromacs_zdist_' + outname + '.pdf', dpi = fig_zdist.dpi)

ax_zpos.set_xlabel('Time (ps)')
ax_zpos.set_ylabel(r'z ($\AA$)')
#ax_zpos.legend(loc = 'lower right')
fig_zpos.tight_layout()
fig_zpos.savefig('gromacs_zpos_' + outname + '.pdf', dpi = fig_zpos.dpi)

ax_helix.set_xlabel('Time (ps)')
ax_helix.set_ylabel('Helical nature (AU)')
#ax_helix.legend(loc = 'lower right')
ax_helix.set_ylim([0.0, 1.1])
fig_helix.tight_layout()
fig_helix.savefig('gromacs_helicity_' + outname + '.pdf', dpi = fig_helix.dpi)

ax_tilt.set_xlabel('Time (ps)')
ax_tilt.set_ylabel('Global tilt (deg)')
#ax_tilt.legend(loc = 'lower right')
ax_tilt.set_ylim([0.0, 180.0])
fig_tilt.tight_layout()
fig_tilt.savefig('gromacs_tilt_' + outname + '.pdf', dpi = fig_tilt.dpi)

