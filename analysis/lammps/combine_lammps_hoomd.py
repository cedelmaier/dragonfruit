#!/usr/bin/env python3

import os
import re
import sys

import freud
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'lib'))
from common import radial_average, ragged_mean
from stylelib.common_styles import septin_runs_stl

def suq_curve(q, N, A, kc, gamma):
    return (N/A) / (kc * q**4 + gamma*q**2)

# Set the style
plt.style.use(septin_runs_stl)

# Read in the LAMMPS CSV file
lammps_csv_filename = "/Users/cedelmaier/Projects/Biophysics/septin_project/supra_cg/dragonfruit/data/20220321/lammps_hoomd_comparison/lammps/nvt_langevin/dumpdata.csv"
df_lammps = pd.read_csv(lammps_csv_filename)
print(df_lammps)

hoomd_csv_filename = "/Users/cedelmaier/Projects/Biophysics/septin_project/supra_cg/dragonfruit/data/20220311/lammps_hoomd_comparison/hoomd/dumpfile_hoomd.csv"
df_hoomd = pd.read_csv(hoomd_csv_filename)
print(df_hoomd)

# Plot everything
x_fft_lammps    = df_lammps['x_fft'].to_numpy()
su_fft_lammps   = df_lammps['su_fft'].to_numpy()
x_fft_lammps    = x_fft_lammps[~np.isnan(x_fft_lammps)]
su_fft_lammps   = su_fft_lammps[~np.isnan(su_fft_lammps)]
lammps_area_mean    = df_lammps['other'].to_numpy()[0]
lammps_nlipids      = df_lammps['other'].to_numpy()[1]

x_fft_hoomd    = df_hoomd['x_fft'].to_numpy()
su_fft_hoomd   = df_hoomd['su_fft'].to_numpy()
x_fft_hoomd    = x_fft_hoomd[~np.isnan(x_fft_hoomd)]
su_fft_hoomd   = su_fft_hoomd[~np.isnan(su_fft_hoomd)]
hoomd_area_mean     = df_hoomd['other'].to_numpy()[0]
hoomd_nlipids       = df_hoomd['other'].to_numpy()[1]

# Plot everything
fig, ax = plt.subplots(1, 1, figsize=(15,10))
# HOOMD and LAMMPS plots
hoomd_scatter   = ax.scatter(x_fft_hoomd[1:], su_fft_hoomd[1:], color = 'b', marker = '+', linewidth = 1, label = 'HOOMD')
lammps_scatter  = ax.scatter(x_fft_lammps[1:], su_fft_lammps[1:], color = 'r', marker = 'o', s = 80, facecolors = 'none', label = 'LAMMPS')
qcutoff = 2.0*np.pi/np.sqrt(hoomd_area_mean)
cutoff_line     = ax.axvline(x = qcutoff, ymin = 0, ymax = 1.0, color = 'k', linestyle = '-', label = 'Frequency cutoff')

# Set up the graphing options
ax.set_ylim(1e0,1e5)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'q ($\sigma^{-1}$)')
ax.set_ylabel(r'$ N \langle | u(q) |^{2} \rangle $ ($\sigma^{2}$)')
ax.legend()

fig.tight_layout()
fig.savefig('lammps_hoomd_comparison_fft_bare.pdf', dpi = fig.dpi)

# Try to do some fits?
from scipy.optimize import curve_fit

# HOOMD fit
idx = np.where(np.greater(x_fft_hoomd, qcutoff))
idx = np.int32(idx[0][0])
jdx = np.where(np.greater(x_fft_hoomd, 1.0))
jdx = np.int32(jdx[0][0])
popt_fft_hoomd_kc, pcov_fft_hoomd_kc    = curve_fit(lambda q, kc: suq_curve(q, hoomd_nlipids, hoomd_area_mean, kc, 0.0), x_fft_hoomd[idx:jdx], su_fft_hoomd[idx:jdx], bounds = ([0.0, np.inf]))

hoomd_fit_kc    = ax.plot(x_fft_hoomd[idx:jdx], suq_curve(x_fft_hoomd[idx:jdx], N = hoomd_nlipids, A = hoomd_area_mean, kc = popt_fft_hoomd_kc[0], gamma = 0.0), color = 'b', linestyle = '--', label = 'HOOMD $k_c$ fit')

# LAMMPS fit
idx = np.where(np.greater(x_fft_lammps, qcutoff))
idx = np.int32(idx[0][0])
jdx = np.where(np.greater(x_fft_lammps, 1.0))
jdx = np.int32(jdx[0][0])
popt_fft_lammps_kc, pcov_fft_lammps_kc  = curve_fit(lambda q, kc: suq_curve(q, lammps_nlipids, lammps_area_mean, kc, 0.0), x_fft_lammps[idx:jdx], su_fft_lammps[idx:jdx], bounds = ([0.0, np.inf]))

lammps_fit_kc   = ax.plot(x_fft_lammps[idx:jdx], suq_curve(x_fft_lammps[idx:jdx], N = lammps_nlipids, A = lammps_area_mean, kc = popt_fft_lammps_kc[0], gamma = 0.0), color = 'r', linestyle = '--', label = 'LAMMPS $k_c$ fit')

print(f"Simulation fits")
print(f"--HOOMD--")
print(f"  kc only:      {popt_fft_hoomd_kc[0]}")
print(f"--LAMMPS--")
print(f"  kc_only:      {popt_fft_lammps_kc[0]}")

ax.legend()

fig.tight_layout()
fig.savefig('lammps_hoomd_comparison_fft_kcfits.pdf', dpi = fig.dpi)
