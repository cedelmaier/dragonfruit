#!/usr/bin/env python3

# XXX: Put a license here

""" Figure 1 hardcoded oneoff script """

import pickle
import os
import subprocess
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
from common import create_datadir

# Set the style for the plots
plt.style.use(septin_poster_stl)

# Move thit into a subclass at some point
def shaded_error(ax, x, y, error, alpha, color, linestyle = 'solid', label = None, linewidth = 1):
    r""" Plot a shaded error bar with colors
    """
    ax.plot(x, y, color = color, linestyle = linestyle, label = label, linewidth = linewidth)
    ax.fill_between(x, y-error, y+error,
                    alpha = alpha, edgecolor = color, facecolor = color, linestyle = linestyle, linewidth = linewidth)

main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data/')
external_path = os.path.abspath('/Volumes/T7/data/septin_project/datasets/')

#simnames = {
#            "neutral_fold_00":      os.path.join(main_path, "rfmonomer_aglipid_11x11_zdepth00_rotx0_50mMKCl_long"),
#            "neutral_fold_15":      os.path.join(external_path, "rfmonomer_aglipid_11x11_zdepth15_50mMKCl"),
#            "neutral_fold_30":      os.path.join(external_path, "rfmonomer_aglipid_11x11_zdepth30_50mMKCl"),
#            "neutral_unfold_15":    os.path.join(external_path, "rfmonomermelt_aglipid_11x11_zdepth15_rotx0_50mMCKl"),
#            "charge_fold_00":       os.path.join(external_path, "gromacs_zdepth00_rotx0_helix_50mMKCl"),
#            "charge_fold_15":       os.path.join(main_path, "agmonomer_aglipid_11x11_zdepth15_rotx0_50mMKCl"),
#            "charge_fold_30":       os.path.join(external_path, "unbiased_50mMKCl_solution"),
#            "charge_unfold_15":     os.path.join(main_path, "agmonomermelt_aglipid_11x11_zdepth15_rotx0_50mMKCl"),
#            }
simnames = {
            "charge_fold_0_rotx0":  os.path.join(main_path, "unbiased_zdepth00_rotx0_helix_50mMKCl"),
            "charge_fold_0_rotx90": os.path.join(main_path, "unbiased_zdepth00_rotx90_helix_50mMKCl"),
           }

seednames = ["N1", "N2", "N3", "N4"]

axis_names = {}
axis_names['zpos']  = r"z ($\AA$)"
axis_names['helix'] = r"Fraction helix (AU)"
axis_names['tilt']  = r"Tilt (deg)"

# Load up the ylow ad yhi dicts
ylow_dict = {"zpos": -30.0,
             "helix": 0.0,
             "tilt": 0.0,
             "pdip": 0.0,
             "pmom": 0.0,
             "hp": 0.0,
             "zforce": -8000.0,
             "perptorque": -800000.0}
yhi_dict = {"zpos": 30.0,
            "helix": 1.05,
            "tilt": 180.0,
            "pdip": 180.0,
            "pmom": 40.0,
            "hp": 180.0,
            "zforce": 8000.0,
            "perptorque": 800000.0}
# Colorblind color cycle
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

# Custom plot options
plt.style.use(septin_paper_stl)

# Colorblind index
cidx = 0

# Just plot the 2 simulations for comparison
simscombined = {
        "flat": os.path.join(simnames['charge_fold_0_rotx0'], 'N1', 'N1.h5'),
        "tilt": os.path.join(simnames['charge_fold_0_rotx90'], 'N2', 'N2.h5'),
        }

fig_zpos, axarr_zpos    = plt.subplots(1, 1, figsize = (4.5, 2.5))
fig_tilt, axarr_tilt    = plt.subplots(1, 1, figsize = (4.5, 2.5))
fig_both, axarr_both    = plt.subplots(2, 1, figsize = (7.5, 5.0), sharex = True)

leaf0_list = []
leaf1_list = []
leaf1_list_top = []
for simname,seedname in simscombined.items():
    print(simname, seedname)

    master_time_df = pd.read_hdf(seedname)
    print(master_time_df)

    # Z position
    z_protein   = master_time_df[['helix_z']].to_numpy().flatten()
    z_leaf0     = master_time_df[['leaflet0_z']].to_numpy().flatten()
    z_leaf1     = master_time_df[['leaflet1_z']].to_numpy().flatten()
    z_lipid     = master_time_df[['lipid_z']].to_numpy().flatten()
    z_pbc       = master_time_df[['unit_cell_z']].to_numpy().flatten()
    # Subtract off the position of the lipid COM from everybody else
    pbc_full_z = z_pbc
    z_protein = z_protein - z_lipid
    z_leaf0 = z_leaf0 - z_lipid
    z_leaf1 = z_leaf1 - z_lipid
    z_pbc = z_pbc/2 - z_lipid
    # Correc the position if under the lower leaflet
    for idx in range(len(z_protein)):
        if z_protein[idx] < z_leaf1[idx]:
            z_protein[idx] = z_pbc[idx] - z_protein[idx]

    xdata = master_time_df.index / 1000.0;

    max_time = xdata[-1]

    # Get the tilt of the helix
    global_tilt = master_time_df[['global_tilt']].to_numpy().flatten()

    # Plot what we want
    plt.figure(fig_zpos)
    ydata = z_protein
    axarr_zpos.plot(xdata, ydata, linewidth = 1, color = CB_color_cycle[cidx])
    axarr_both[0].plot(xdata, ydata, linewidth = 1, color = CB_color_cycle[cidx])

    leaf0_list.append(z_leaf0)
    leaf1_list.append(z_leaf1)
    leaf1_list_top.append(z_leaf1 + pbc_full_z)

    plt.figure(fig_tilt)
    ydata = global_tilt
    axarr_tilt.plot(xdata, ydata, linewidth = 1, color = CB_color_cycle[cidx])
    axarr_both[1].plot(xdata, ydata, linewidth = 1, color = CB_color_cycle[cidx])

    cidx += 1

# Do the leaflet limits
plt.figure(fig_zpos)
leaf0_mean = np.mean(np.array(leaf0_list), axis=0)
leaf1_mean = np.mean(np.array(leaf1_list), axis=0)
leaf1_top_mean = np.mean(np.array(leaf1_list_top), axis=0)
leaf0_std = np.std(np.array(leaf0_list), axis=0, ddof=1)
leaf1_std = np.std(np.array(leaf1_list), axis=0, ddof=1)
leaf1_top_std = np.std(np.array(leaf1_list_top), axis=0, ddof=1)

# Shaded error bars for this
shaded_error(axarr_zpos, xdata, leaf0_mean, leaf0_std, alpha = 0.5, color = 'slategrey')
shaded_error(axarr_zpos, xdata, leaf1_mean, leaf1_std, alpha = 0.5, color = 'slategrey')
shaded_error(axarr_zpos, xdata, leaf1_top_mean, leaf1_top_std, alpha = 0.5, color = 'slategrey')
# Alsof for the combined version
shaded_error(axarr_both[0], xdata, leaf0_mean, leaf0_std, alpha = 0.5, color = 'slategrey')
shaded_error(axarr_both[0], xdata, leaf1_mean, leaf1_std, alpha = 0.5, color = 'slategrey')
shaded_error(axarr_both[0], xdata, leaf1_top_mean, leaf1_top_std, alpha = 0.5, color = 'slategrey')

# Set the limits
axarr_zpos.set_ylim(ylow_dict['zpos'], yhi_dict['zpos'])
axarr_zpos.set_xlim(0.0, 200.0)

# Labels and such
axarr_zpos.set_xlabel(r"Time (ns)")
axarr_zpos.set_ylabel(axis_names['zpos'])

plt.figure(fig_zpos)
fig_zpos.tight_layout()
plt.savefig("flat_tilt_zpos_comparison.png", dpi = 600)

# Now do the tilt
plt.figure(fig_tilt)
axarr_tilt.set_ylim(ylow_dict['tilt'], yhi_dict['tilt'])
axarr_tilt.set_xlim(0.0, 200.0)

axarr_tilt.set_xlabel(r"Time(ns)")
axarr_tilt.set_ylabel(axis_names['tilt'])

axarr_tilt.legend(["Flat", "Tilted"], loc = 'lower right')

fig_tilt.tight_layout()
plt.savefig("flat_tilt_tilt_comparison.png", dpi = 600)

# Do the combination
plt.figure(fig_both)
axarr_both[1].set_ylim(ylow_dict['tilt'], yhi_dict['tilt'])
axarr_both[0].set_ylim(ylow_dict['zpos'], yhi_dict['zpos'])

axarr_both[1].set_xlabel(r"Time(ns)")
axarr_both[1].set_ylabel(axis_names['tilt'])
axarr_both[0].set_ylabel(axis_names['zpos'])

axarr_both[1].legend(["Flat", "Tilted"], loc = 'lower right')

fig_both.tight_layout()
plt.savefig("flat_tilt_both.png", dpi = 600)
