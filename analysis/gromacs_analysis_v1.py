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

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src', 'lib'))
from stylelib.common_styles import septin_runs_stl

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
                    'coiled/zdepth_00angstroms/s1',
                    'coiled/zdepth_10angstroms/s1',
                    'coiled/zdepth_20angstroms/s1',
                    'coiled/zdepth_30angstroms/s1',
                    #'unfolded/zdepth_00angstroms/s1',
                    #'unfolded/zdepth_10angstroms/s1',
                    #'unfolded/zdepth_10angstroms/s2',
                    #'unfolded/zdepth_10angstroms/s3',
                    #'unfolded/zdepth_10angstroms/s4',
                    #'unfolded/zdepth_20angstroms/s1',
                    #'unfolded/zdepth_20angstroms/s2',
                    #'unfolded/zdepth_20angstroms/s3',
                    #'unfolded/zdepth_20angstroms/s4',
                    #'unfolded/zdepth_30angstroms/s1',
                    ]

# Set up an output name
outname = 'coiled_all'
#outname = 'unfolded_all'
#outname = 'unfolded_10'
#outname = 'unfolded_20'

# Set up the plots beforehand
plt.style.use(septin_runs_stl)
fig_zdist, ax_zdist = plt.subplots(1, 1, figsize = (15, 10))
fig_twist, ax_twist = plt.subplots(1, 1, figsize = (15, 10))
fig_resid, ax_resid = plt.subplots(1, 1, figsize = (15, 10))

# Set up a stride and only graph every stride points
stride = 10

# For each simulation, loop through and get the relevant data and plot it
for sname in simulation_names:
    filepath = os.path.join(main_path, sname)
    filenames = filepath.split('/')[-1]

    hd5_filename = os.path.join(filepath, filenames + '.h5')
    pkl_filename = os.path.join(filepath, filenames + '.pickle')

    master_df = pd.read_hdf(hd5_filename)
    helix_analysis = None
    with open(pkl_filename, 'rb') as f:
        helix_analysis = pickle.load(f)

    # Get the times for this
    times = master_df.index

    # Generate the depth and name for the label
    coilname = sname.split('/')[0]
    depth = sname.split('_')[1][0:2]
    seedname = sname.split('/')[-1]
    labelname = coilname + ' +' + depth + ' ' + seedname

    # Compute zdist plot
    zdist = np.abs(master_df['helix_z'] - master_df['lipid_z'])
    ax_zdist.plot(times, zdist, label = labelname)

    # Get the average twist
    avg_twist = helix_analysis.results.local_twists.mean(axis = 1)
    ax_twist.plot(times, avg_twist, label = labelname)

    # Get the residues per turn
    nres_per_turn = helix_analysis.results.local_nres_per_turn.mean(axis = 1)
    ax_resid.plot(times, nres_per_turn, label = labelname)


ax_zdist.set_xlabel('Time (ps)')
ax_zdist.set_ylabel('Z distance (Angstroms)')
ax_zdist.legend(loc = 'lower right')
#ax_zdist.set_ylim([0.0, 60.0])
fig_zdist.tight_layout()
fig_zdist.savefig('gromacs_zdist_' + outname + '.pdf', dpi = fig_zdist.dpi)

ax_twist.set_xlabel('Time (ps)')
ax_twist.set_ylabel('Average twist (degrees)')
ax_twist.legend()
fig_twist.tight_layout()
fig_twist.savefig('gromacs_twist_' + outname + '.pdf')

ax_resid.set_xlabel('Time (ps)')
ax_resid.set_ylabel('Average residues per turn')
ax_resid.legend()
fig_resid.tight_layout()
fig_resid.savefig('gromacs_nres_'+ outname + '.pdf')
