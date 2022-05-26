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

main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/')
#simulation_names = [
#                    #'ahcoil_00angstrom_v1',
#                    #'ahcoil_10angstrom_v1',
#                    #'ahcoil_20angstrom_v1',
#                    #'ahcoil_30angstrom_v1',
#                    'ahhelix_00angstrom_v1',
#                    'ahhelix_10angstrom_v1',
#                    'ahhelix_20angstrom_v1',
#                    #'ahhelix_30angstrom_v1',
#                    ]
simulation_names = [
                    #'coiled/zdepth_00angstroms/s1',
                    #'coiled/zdepth_10angstroms/s1',
                    #'coiled/zdepth_20angstroms/s1',
                    #'coiled/zdepth_30angstroms/s1',
                    'unfolded/zdepth_00angstroms/s1',
                    'unfolded/zdepth_10angstroms/s1',
                    #'unfolded/zdepth_10angstroms/s2',
                    #'unfolded/zdepth_10angstroms/s3',
                    #'unfolded/zdepth_10angstroms/s4',
                    'unfolded/zdepth_20angstroms/s1',
                    #'unfolded/zdepth_20angstroms/s2',
                    #'unfolded/zdepth_20angstroms/s3',
                    #'unfolded/zdepth_20angstroms/s4',
                    'unfolded/zdepth_30angstroms/s1',
                    ]
simulation_legends = [
                      #r'Helix, 0$\AA$ depth',
                      #r'Helix, 10$\AA$ depth',
                      #r'Helix, 20$\AA$ depth',
                      #r'Helix, 30$\AA$ depth',
                      r'Unfolded, 0$\AA$ depth',
                      r'Unfolded, 10$\AA$ depth, s1',
                      #r'Unfolded, 10$\AA$ depth, s2',
                      #r'Unfolded, 10$\AA$ depth, s3',
                      #r'Unfolded, 10$\AA$ depth, s4',
                      r'Unfolded, 20$\AA$ depth, s1',
                      #r'Unfolded, 20$\AA$ depth, s2',
                      #r'Unfolded, 20$\AA$ depth, s3',
                      #r'Unfolded, 20$\AA$ depth, s4',
                      r'Unfolded, 30$\AA$ depth',
                     ]


# Set up an output name
#outname = 'coiled_all'
outname = 'unfolded_main'
#outname = 'unfolded_10'
#outname = 'unfolded_20'

# Set up the plots beforehand
plt.style.use(septin_poster_stl)
fig_zdist, ax_zdist = plt.subplots(1, 1, figsize = (15, 10))
fig_zpos , ax_zpos  = plt.subplots(1, 1, figsize = (15, 10))
fig_twist, ax_twist = plt.subplots(1, 1, figsize = (15, 10))
fig_resid, ax_resid = plt.subplots(1, 1, figsize = (15, 10))
fig_helix, ax_helix = plt.subplots(1, 1, figsize = (15, 10))

# Set up a stride and only graph every stride points
stride = 10

# For each simulation, loop through and get the relevant data and plot it
cidx = 4
for idx,sname in enumerate(simulation_names):
    filepath = os.path.join(main_path, sname)
    filenames = filepath.split('/')[-1]

    hd5_filename = os.path.join(filepath, filenames + '.h5')
    #pkl_filename = os.path.join(filepath, filenames + '.pickle')

    master_df = pd.read_hdf(hd5_filename)
    #helix_analysis = None
    #with open(pkl_filename, 'rb') as f:
    #    helix_analysis = pickle.load(f)

    # Get the times for this
    times = master_df.index

    # Generate the depth and name for the label
    coilname = sname.split('/')[0]
    depth = sname.split('_')[1][0:2]
    seedname = sname.split('/')[-1]
    #labelname = coilname + ' +' + depth + ' ' + seedname
    labelname = simulation_legends[idx]

    # Compute zdist plot
    #zdist = np.abs(master_df['helix_z'] - master_df['lipid_z'])
    #ax_zdist.plot(times, zdist, label = labelname)

    ## Get the average twist
    #avg_twist = helix_analysis.results.local_twists.mean(axis = 1)
    #ax_twist.plot(times, avg_twist, label = labelname)

    ## Get the residues per turn
    #nres_per_turn = helix_analysis.results.local_nres_per_turn.mean(axis = 1)
    #ax_resid.plot(times, nres_per_turn, label = labelname)

    # Compute the z position
    z_pos = master_df[['helix_z']].to_numpy().flatten()
    leaflet0_pos = master_df[['leaflet0_z']].to_numpy().flatten()
    leaflet1_pos = master_df[['leaflet1_z']].to_numpy().flatten()
    lipid_pos = master_df[['lipid_z']].to_numpy().flatten()
    z_pbc = master_df[['unit_cell_z']].to_numpy().flatten()
    # Subtract off the lipid COM position
    z_pos = z_pos - lipid_pos
    leaflet0_pos = leaflet0_pos - lipid_pos
    leaflet1_pos = leaflet1_pos - lipid_pos
    z_pbc = z_pbc/2 - lipid_pos
    # Correct the position if under the lower leaflet
    for idx in range(len(z_pos)):
        if z_pos[idx] < leaflet1_pos[idx]:
            z_pos[idx] = z_pbc[idx] - z_pos[idx]

    # Smooth these in the same way
    z_pos_hat = savgol_filter(z_pos, 51, 3)
    leaflet0_hat = savgol_filter(leaflet0_pos, 51, 3)
    leaflet1_hat = savgol_filter(leaflet1_pos, 31, 3)
    #ax_zpos.plot(times[::stride], z_pos[::stride], label = labelname)
    #ax_zpos.plot(times[::stride], leaflet0_pos[::stride], color = 'k')
    #ax_zpos.plot(times[::stride], leaflet1_pos[::stride], color = 'k')
    ax_zpos.plot(times, z_pos_hat, label = labelname, color = CB_color_cycle[cidx])
    ax_zpos.plot(times, leaflet0_hat, color = 'k')
    ax_zpos.plot(times, leaflet1_hat, color = 'k')

    # Compute the helicity
    helicity = master_df[['helicity']].to_numpy().flatten()
    # Smooth this data, as it is gross right now
    yhat = savgol_filter(helicity, 51, 3)
    yhat[yhat <= 0.0] = 0.0
    yhat[yhat >= 1.0] = 1.0
    ax_helix.plot(times, yhat, label = labelname, color = CB_color_cycle[cidx])

    cidx += 1


#ax_zdist.set_xlabel('Time (ps)')
#ax_zdist.set_ylabel('Z distance (Angstroms)')
#ax_zdist.legend(loc = 'lower right')
##ax_zdist.set_ylim([0.0, 60.0])
#fig_zdist.tight_layout()
#fig_zdist.savefig('gromacs_zdist_' + outname + '.pdf', dpi = fig_zdist.dpi)

#ax_twist.set_xlabel('Time (ps)')
#ax_twist.set_ylabel('Average twist (degrees)')
#ax_twist.legend()
#fig_twist.tight_layout()
#fig_twist.savefig('gromacs_twist_' + outname + '.pdf')
#
#ax_resid.set_xlabel('Time (ps)')
#ax_resid.set_ylabel('Average residues per turn')
#ax_resid.legend()
#fig_resid.tight_layout()
#fig_resid.savefig('gromacs_nres_'+ outname + '.pdf')

ax_zpos.set_xlabel('Time (ps)')
ax_zpos.set_ylabel(r'z ($\AA$)')
#ax_zpos.legend()
ax_zpos.set_ylim([-25.0, 75.0])
fig_zpos.tight_layout()
fig_zpos.savefig('gromacs_zpos_' + outname + '.pdf')

ax_helix.set_xlabel('Time (ps)')
ax_helix.set_ylabel('Helical nature (AU)')
#ax_helix.legend()
ax_helix.set_ylim([0.0, 1.1])
fig_helix.tight_layout()
fig_helix.savefig('gromacs_helicity_' + outname + '.pdf')
