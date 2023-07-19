#!/usr/bin/env python3

# XXX: Put a license here

""" Figure 3 hardcoded oneoff script """

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
import plumed
from scipy.signal import savgol_filter

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src', 'lib'))
from stylelib.common_styles import *
from common import create_datadir

# Set the style for the plots
plt.style.use(septin_poster_stl)

def adjust_fes(dz, alpha, fes, start_index):
    r""" Adjust the FES and everything based on the start index to slice off some rows in alpha
    """
    dz_adj       = dz[start_index:][:]
    alpha_adj    = alpha[start_index:][:]
    fes_adj      = fes[start_index:][:]

    mmin = np.min(fes_adj)

    return dz_adj, alpha_adj, fes_adj, mmin

def plot_fes(dz, alpha, fes, fname):
    r""" Plot the FES with the given name (x and y axis will always be the same
    """
    fig, ax = plt.subplots(1, 1, figsize = (10, 10))
    ax.contour(dz, alpha, fes, levels=range(0,100,10), linewidths=0.5, colors='k')
    cntr = ax.contourf(dz, alpha, fes, levels=range(0,100), cmap='rainbow')
    plt.colorbar(cntr, label="FES [kJ/mol]")
    
    ax.set_xlim(0.0, 5.0)
    ax.set_ylim(0.0, 13.1)
    
    ax.set_xlabel(r'z (nm)')
    ax.set_ylabel(r'Helical content (AU)')
    
    fig.tight_layout()
    plt.savefig(fname, dpi = fig.dpi)
    plt.close()

# Read in the free energy landscape from the sum_hills utility
main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data/enhanced_sampling/metadynamics/unbias_neutral')

data_fes_original = plumed.read_as_pandas(os.path.join(main_path, 'fes_sumhills_neutral_entirespace.dat'))
npoints = 256
# The total size of the run for later
totalrunsize = 1455128

dz_original    = np.array(data_fes_original["dz"]).reshape(npoints, npoints)
alpha_original = np.array(data_fes_original["alpha"]).reshape(npoints, npoints)
fes_original   = np.array(data_fes_original["file.free"]).reshape(npoints, npoints)

fig, ax = plt.subplots(1, 1, figsize = (10, 10))

# Slice off the first two rows of alpha (0 and whatever comes after) to get close to 0.1
# This also involves slicing off Z and fes_original
start_index = 4
dz_adjusted, alpha_adjusted, fes_adjusted, min_adjusted = adjust_fes(dz_original, alpha_original, fes_original, start_index)

# Readjust the minimum to zero based on the min of fes_adjusted
print(f"New minimum in adjusted[{start_index}] space, {min_adjusted}")
fes_adjusted = fes_adjusted - min_adjusted

plot_fes(dz_adjusted, alpha_adjusted, fes_adjusted, "neutral_fes_adjusted.pdf")

# XXX TODO add the X1-4 points on the graph directly

########
# Convergence testing
# look at ~/plot_convergence scripts for original information
########

# Biased energy landscape
data_biased = plumed.read_as_pandas(os.path.join(main_path, 'ff_z_alpha_biased.dat'))
dz_biased       = np.array(data_biased["dz"]).reshape(npoints, npoints)
alpha_biased    = np.array(data_biased["alpha"]).reshape(npoints, npoints)
fes_biased      = np.array(data_biased["ff_z_alpha_biased"]).reshape(npoints, npoints)
# Adjust the same as the original
dz_biased_adj, alpha_biased_adj, fes_biased_adj, min_biased_adjusted = adjust_fes(dz_biased, alpha_biased, fes_biased, start_index)
fes_biased_adj = fes_biased_adj - min_biased_adjusted
print(f"New minimum in biased adjusted[{start_index}] space, {min_biased_adjusted}")

plot_fes(dz_biased_adj, alpha_biased_adj, fes_biased_adj, "neutral_fes_biasedadjusted.pdf")

# Normal reweighting scheme (exponential)
data_as = plumed.read_as_pandas(os.path.join(main_path, "ff_z_alpha_as.dat"))
dz_biased_as    = np.array(data_as["dz"]).reshape(npoints, npoints)
alpha_biased_as = np.array(data_as["alpha"]).reshape(npoints, npoints)
fes_biased_as   = np.array(data_as["ff_z_alpha_as"]).reshape(npoints, npoints)
# Adjust same
dz_biased_as_adj, alpha_biased_as_adj, fes_biased_as_adj, min_biased_as_adj = adjust_fes(dz_biased_as, alpha_biased_as, fes_biased_as, start_index)
fes_biased_as_adj = fes_biased_as_adj - min_biased_as_adj
print(f"New minimum in reweight biased adjusted[{start_index}] space, {min_biased_as_adj}")

plot_fes(dz_biased_as_adj, alpha_biased_as_adj, fes_biased_as_adj, "neutral_fes_as_biasedadjusted.pdf")

# Tiwary reweighting
data_tw = plumed.read_as_pandas(os.path.join(main_path, "ff_z_alpha_tw.dat"))
dz_biased_tw    = np.array(data_tw["dz"]).reshape(npoints, npoints)
alpha_biased_tw = np.array(data_tw["alpha"]).reshape(npoints, npoints)
fes_biased_tw   = np.array(data_tw["ff_z_alpha_tw"]).reshape(npoints, npoints)
# Adjust same
dz_biased_tw_adj, alpha_biased_tw_adj, fes_biased_tw_adj, min_biased_tw_adj = adjust_fes(dz_biased_tw, alpha_biased_tw, fes_biased_tw, start_index)
fes_biased_tw_adj = fes_biased_tw_adj - min_biased_tw_adj
print(f"New minimum in reweight (Tiwary) biased adjusted[{start_index}] space, {min_biased_tw_adj}")

plot_fes(dz_biased_tw_adj, alpha_biased_tw_adj, fes_biased_tw_adj, "neutral_fes_tw_biasedadjusted.pdf")

