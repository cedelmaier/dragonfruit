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

def setup_broken_axis(axarr):
    r""" Setup a broken plot for ourselves
    """
    axarr[0].spines['right'].set_visible(False)
    axarr[0].tick_params(labelright = False)
    axarr[1].spines['left'].set_visible(False)
    axarr[1].tick_params(labelleft = False)
    axarr[1].yaxis.tick_right()
    axarr[1].tick_params(right=False)

    # Create the breaks
    d = 0.01

    # Gross
    kwargs = dict(transform=axarr[0].transAxes, color = 'k', clip_on = False)
    axarr[0].plot((1-d, 1+d), (-d, d), linewidth=1, **kwargs)
    axarr[0].plot((1-d, 1+d), (1-d, 1+d), linewidth=1, **kwargs)
    kwargs.update(transform=axarr[1].transAxes)
    axarr[1].plot((-d, d), (-d, d), linewidth=1, **kwargs)
    axarr[1].plot((-d, d), (1-d, 1+d), linewidth=1, **kwargs)

main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data/figures/figure1')

# Create some external data directories
combined_datadir = create_datadir(os.getcwd(), datadir_name = "combined")
combinedsims = {
        "fold_00": ["neutral_fold_00", "charge_fold_00"],
        "fold_15": ["neutral_fold_15", "charge_fold_15"],
        "fold_30": ["neutral_fold_30", "charge_fold_30"],
        "unfold_15": ["neutral_unfold_15", "charge_unfold_15"],
        }
chargetypes = ["neutral", "charge"]
membrane_colors = ["slategrey", "silver"]

axis_names = {}
#axis_names['zpos']  = r"z ($\AA$)"
#axis_names['helix'] = r"Fraction helix (AU)"
axis_names['zpos'] = r"z (nm)"
axis_names['alpharmsd'] = r"Helical content (AU)"

# Load up the ylow ad yhi dicts
ylow_dict = {"zpos": 0.0,
             "alpharmsd": 0.0,
             "helix": 0.0,
             "tilt": 0.0,
             "pdip": 0.0,
             "pmom": 0.0,
             "hp": 0.0,
             "zforce": -8000.0,
             "perptorque": -800000.0}
yhi_dict = {"zpos": 10.0,
            "alpharmsd": 13.1,
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

for simname,combinedfiles in combinedsims.items():
    print(simname, combinedfiles)

    # Create the figures
    fig_zpos_mean, axarr_zpos_mean  = plt.subplots(1, 2, figsize = (5.0, 2.5), sharey = True)
    fig_heli_mean, axarr_heli_mean  = plt.subplots(1, 2, figsize = (5.0, 2.5), sharey = True)

    # Load the HDF5 files
    fname_neutral   = combinedfiles[0] + ".h5"
    fname_charge    = combinedfiles[1] + ".h5"

    dfs = {}
    for chargetype in chargetypes:
        if chargetype == "neutral":
            dfs[chargetype] = pd.read_hdf(fname_neutral)
        else:
            dfs[chargetype] = pd.read_hdf(fname_charge)

    # This is very gross and hardcoded
    for itype in range(2):
        ctype = chargetypes[itype]
        membrane_color = membrane_colors[itype]
        traj_color = CB_color_cycle[itype]

        xdata           = dfs[ctype]['timepoints'].to_numpy().flatten()

        max_time = xdata[-1]

        leaf0_mean      = dfs[ctype]['leaf0_mean'].to_numpy().flatten()
        leaf1_mean      = dfs[ctype]['leaf1_mean'].to_numpy().flatten()
        leaf1_top_mean  = dfs[ctype]['leaf1_top_mean'].to_numpy().flatten()
        center_mean     = dfs[ctype]['center_mean'].to_numpy().flatten()
        heli_mean       = dfs[ctype]['heli_mean'].to_numpy().flatten()

        leaf0_std       = dfs[ctype]['leaf0_std'].to_numpy().flatten()
        leaf1_std       = dfs[ctype]['leaf1_std'].to_numpy().flatten()
        leaf1_top_std   = dfs[ctype]['leaf1_top_std'].to_numpy().flatten()
        center_std      = dfs[ctype]['center_std'].to_numpy().flatten()
        heli_std        = dfs[ctype]['heli_std'].to_numpy().flatten()

        # Z position
        # Draw the membranes
        mlabel = None
        for iax in range(2):
            # Only plot membranes from the neutralversion 
            if ctype == "charge":
                continue
            if simname == "unfold_15":
                mlabel = "Membrane"
            else:
                mlabel = None
            shaded_error(axarr_zpos_mean[iax], xdata, leaf0_mean, leaf0_std, alpha = 0.5, color = membrane_color)
            shaded_error(axarr_zpos_mean[iax], xdata, leaf1_mean, leaf1_std, alpha = 0.5, color = membrane_color, label = mlabel)
            shaded_error(axarr_zpos_mean[iax], xdata, leaf1_top_mean, leaf1_top_std, alpha = 0.5, color = membrane_color)

        # Draw the actual trajectory
        for iax in range(2):
            print(f"{simname}, {iax}")
            if (simname == "unfold_15" and ctype == "neutral"):
                mlabel = "Neutral cap"
            elif (simname == "unfold_15" and ctype == "charge"):
                mlabel = "Charge cap"
            else:
                mlabel = None
            shaded_error(axarr_zpos_mean[iax], xdata, center_mean, center_std, alpha = 0.5, color = traj_color, label = mlabel)

        # Helices
        # This is much simpler as we don't have other floppy bits
        for iax in range(2):
            shaded_error(axarr_heli_mean[iax], xdata, heli_mean, heli_std, alpha = 0.5, color = traj_color)

    ########
    # Z position
    ########
    # Swap to the proper figure handle
    plt.figure(fig_zpos_mean)
    # Set the limits
    axarr_zpos_mean[0].set_ylim(ylow_dict["zpos"], yhi_dict["zpos"])
    axarr_zpos_mean[1].set_ylim(ylow_dict["zpos"], yhi_dict["zpos"])
    axarr_zpos_mean[0].set_xlim(0.0, 10.0)
    axarr_zpos_mean[1].set_xlim(max_time - 100.0, max_time)

    # Set the plot breaks
    setup_broken_axis(axarr_zpos_mean)

    # Create the Z-labels
    spaceshift = 0.10
    fig_zpos_mean.supxlabel(r"Time (ns)", x = 0.5+spaceshift/2.0)
    axarr_zpos_mean[0].set_ylabel(axis_names["zpos"])

    # Output the Z-axis
    plt.figure(fig_zpos_mean)
    if simname == "unfold_15":
        axarr_zpos_mean[0].legend()
    fig_zpos_mean.tight_layout()
    plt.subplots_adjust(wspace=spaceshift, bottom=0.16)
    output_filename = combined_datadir + "/z_" + simname + ".png"
    plt.savefig(output_filename, dpi = 600)
    plt.close()

    ########
    # Helical content
    ########
    plt.figure(fig_heli_mean)
    axarr_heli_mean[0].set_ylim(ylow_dict["alpharmsd"], yhi_dict["alpharmsd"])
    axarr_heli_mean[1].set_ylim(ylow_dict["alpharmsd"], yhi_dict["alpharmsd"])
    axarr_heli_mean[0].set_xlim(0.0, 10.0)
    axarr_heli_mean[1].set_xlim(max_time - 100.0, max_time)

    # Set the plot breaks
    setup_broken_axis(axarr_heli_mean)

    # Create the labels
    spaceshift = 0.10
    fig_heli_mean.supxlabel(r"Time (ns)", x = 0.5+spaceshift/2.0)
    axarr_heli_mean[0].set_ylabel(axis_names["alpharmsd"])

    # Output helical conctent
    plt.figure(fig_heli_mean)
    fig_heli_mean.tight_layout()
    plt.subplots_adjust(wspace=spaceshift, bottom=0.16)
    output_filename = combined_datadir + "/heli_" + simname + ".png"
    plt.savefig(output_filename, dpi = 600)
    plt.close()

