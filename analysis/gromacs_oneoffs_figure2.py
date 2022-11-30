#!/usr/bin/env python3

# XXX: Put a license here

""" Figure 2 hardcoded oneoff script """

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

import plumed

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src', 'lib'))
from stylelib.common_styles import *
from dragonfruit_io import read_xvg

# Move thit into a subclass at some point
def shaded_error(ax, x, y, error, alpha, color, linestyle = 'solid', label = None, linewidth = 1):
    r""" Plot a shaded error bar with colors
    """
    ax.plot(x, y, color = color, linestyle = linestyle, label = label, linewidth = linewidth)
    ax.fill_between(x, y-error, y+error,
                    alpha = alpha, edgecolor = color, facecolor = color, linestyle = linestyle, linewidth = linewidth)

main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data/enhanced_sampling')

umbrella_names = {
        "umbrella_outsidein_200kJ": "rfmonomer_aglipid_11x11_zdepth80_50mMKCl_smd",
        "umbrella_outsidein_1000kJ": "rfmonomer_aglipid_11x11_zdepth80_50mMKCl_smd_higherkappa"
        }

umbrella_labels = {
        "umbrella_outsidein_200kJ": r"200 kJ mol$^{-1}$ nm$^{-2}$ SMD",
        "umbrella_outsidein_1000kJ": r"1000 kJ mol$^{-1}$ nm$^{-2}$ SMD",
        }

metadynamics_names = {
        "metad_20k": "rfmonomer_aglipid_11x11_zdepth00_50mMKCl_metad",
        "metad_40k": "rfmonomer_aglipid_11x11_zdepth00_50mMKCl_metad_v2"
        }

# Colorblind color cycle
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

# Custom plot options
plt.style.use(septin_paper_stl)

dfs = {}

# Read in the XVG file for the umbrella sampling. Note, there are several of these
cidx = 0
fig_umbrella_all, ax_umbrella_all = plt.subplots(1, 1, figsize=(5.0,2.50))
for simname,seedname in umbrella_names.items():
    print(simname, seedname)

    # Look for an XVG file
    profile_xvg_filename = os.path.join(main_path, seedname, "profile.xvg")
    [metadata, num_data, df] = read_xvg(profile_xvg_filename)

    fig_umbrella, ax_umbrella = plt.subplots(1, 1, figsize=(5.0, 2.5))

    # Plot the free energy
    df.columns = ['FE']
    df.index.name = "z"
    xdata = df.index.to_numpy().flatten()*10.0
    ydata = df[['FE']].to_numpy().flatten()
    ax_umbrella.plot(xdata, ydata, linewidth = 1, color = CB_color_cycle[cidx])

    # Set limits
    ax_umbrella.set_xlim(0.0, 100.0)

    ax_umbrella.set_xlabel(r"z ($\AA$)")
    ax_umbrella.set_ylabel(r"Free energy (kJ mol$^{-1}$)")
    
    fig_umbrella.tight_layout()
    fig_umbrella.savefig("figure2_umbrella_" + simname + ".pdf", dpi = fig_umbrella.dpi)

    # Add to our overall dataframe mess
    dfs[simname] = df
    plt.close()

    ydata = ydata - np.min(ydata)
    ax_umbrella_all.plot(xdata, ydata, linewidth = 1, color = CB_color_cycle[cidx], label = umbrella_labels[simname])

    cidx += 1

plt.figure(fig_umbrella_all)
ax_umbrella_all.set_xlim(0.0, 100.0)
ax_umbrella_all.set_xlabel(r"z ($\AA$)")
ax_umbrella_all.set_ylabel(r"Free energy (kJ mol$^{-1}$)")
ax_umbrella_all.legend()
fig_umbrella_all.tight_layout()
fig_umbrella_all.savefig("figure2_umbrella_comparison.pdf", dpi = fig_umbrella_all.dpi)
plt.close()

cidx = 0
fig_metad_all, ax_metad_all = plt.subplots(1, 1, figsize=(5.0,2.5))
# What about metadynamics?
for simname,seedname in metadynamics_names.items():
    print(simname, seedname)

    # Load the colvar file
    df_colvar = plumed.read_as_pandas(os.path.join(main_path, seedname, "colvar_metad.dat"))

    # Plot the position over time, and alpha vs. position
    fig_position, ax_position = plt.subplots(1, 1, figsize=(5.0,2.5))
    ax_position.plot(df_colvar["time"], df_colvar["z_dist.z"], 'o', ms=1)
    ax_position.set_xlabel(r"Simulation frame")
    ax_position.set_ylabel(r"z ($\AA)")
    fig_position.tight_layout()
    fig_position.savefig("figure2_metad_" + simname + "_zpos.pdf", dpi = fig_position.dpi)
    plt.close()

    # Position vs alpha
    fig_z_alpha, ax_z_alpha = plt.subplots(1, 1, figsize=(5.0,2.5))
    ax_z_alpha.plot(df_colvar["z_dist.z"], df_colvar["alpha"], 'o', ms=2, label="MetaD")
    ax_z_alpha.set_xlabel(r"z ($\AA$)")
    ax_z_alpha.set_ylabel(r"Helical residues")
    fig_z_alpha.tight_layout()
    fig_z_alpha.savefig("figure2_metad_" + simname + "_zvsalpha.pdf", dpi = fig_z_alpha.dpi)
    plt.close()

    df_fes = plumed.read_as_pandas(os.path.join(main_path, seedname, "fes_100.dat"))
    if simname == "metad_40k":
        df_fes = plumed.read_as_pandas(os.path.join(main_path, seedname, "fes_200.dat"))
    fig_fes, ax_fes = plt.subplots(1, 1, figsize=(5.0,2.5))
    ax_fes.plot(df_fes["z_dist.z"]*10.0, df_fes["file.free"], linewidth = 1, color = CB_color_cycle[cidx])
    ax_fes.set_xlabel(r"z ($\AA$)")
    ax_fes.set_ylabel(r"Free Energy (kJ mol$^{-1}$)")
    fig_fes.tight_layout()
    fig_fes.savefig("figure2_metad_" + simname + "_fes.pdf", dpi = fig_fes.dpi)
    plt.close()

    ax_metad_all.plot(df_fes["z_dist.z"]*10.0, df_fes["file.free"], linewidth = 1, color = CB_color_cycle[cidx], label = simname)

    df_new = df_fes[['z_dist.z', 'file.free']].copy(deep=True)
    df_new.columns = ['z', 'FE']
    df_new.set_index('z', inplace=True, drop=True)

    dfs[simname] = df_new

    cidx += 1

plt.figure(fig_metad_all)
ax_metad_all.set_xlim(0.0, 100.0)
ax_metad_all.set_xlabel(r"z ($\AA$)")
ax_metad_all.set_ylabel(r"Free energy (kJ mol$^{-1}$)")
ax_metad_all.legend()
fig_metad_all.tight_layout()
fig_metad_all.savefig("figure2_metad_comparison.pdf", dpi = fig_metad_all.dpi)
plt.close()

# Now plot the different energies!
fig_fes, ax_fes = plt.subplots(1, 1, figsize = (5.0, 2.5))
cidx = 0
for dfname,df in dfs.items():
    # Convert nm to angstroms
    xdata = df.index.to_numpy().flatten()*10.0

    # Do the ydata as well, find minimum shift, etc
    ydata = df[['FE']].to_numpy().flatten()
    # Find the minimum and shift to zero
    ydata = ydata - np.min(ydata)

    ax_fes.plot(xdata, ydata, linewidth=1, color = CB_color_cycle[cidx], label = dfname)

    cidx += 1

# Now set limits, etc
ax_fes.set_xlim(-50.0, 100.0)

ax_fes.set_xlabel(r"z ($\AA$)")
ax_fes.set_ylabel(r"Free energy (kJ mol$^{-1}$)")

fig_fes.legend()

fig_fes.tight_layout()
fig_fes.savefig("figure2_fes_comparison.pdf", dpi = fig_fes.dpi)
plt.close()


