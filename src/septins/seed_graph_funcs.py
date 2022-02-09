# Something

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from common import moving_average

def graph_seed_scatter(sdlabel,
                       xdata,
                       ydata,
                       ax,
                       mtitle = '',
                       xtitle = '',
                       ytitle = '',
                       color = 'b',
                       xlabel = True):
                  
    r""" Generic 1D plotting function for data
    """
    ax.set_title(mtitle)
    ax.set_ylabel(ytitle)
    if xlabel: ax.set_xlabel(xtitle)
    ax.scatter(xdata, ydata, color = color, label = sdlabel)

def graph_seed_lifetime_distribution(sdlabel,
                                     df,
                                     axarr,
                                     mtitle = '',
                                     xtitle = '',
                                     ytitle = '',
                                     color = 'b',
                                     xlabel = True):

    r""" Plotting function for total lifetime of septin attachments
    """
    hist_free, bin_edges_free = np.histogram(df['free'], bins = 10, range = (0, 200))
    hist_near, bin_edges_near = np.histogram(df['near'], bins = 10, range = (0, 200))
    hist_surf, bin_edges = np.histogram(df['surface'], bins = 10, range = (0, 200))
    hist_inte, bin_edges = np.histogram(df['intermediate'], bins = 10, range = (0, 200))
    hist_deep, bin_edges = np.histogram(df['deep'], bins = 10, range = (0, 200))
    bin_mids = moving_average(bin_edges)
    bin_width = bin_mids[1] - bin_mids[0]

    axarr[0][0].bar(bin_mids, hist_free, bin_width)
    axarr[0][1].bar(bin_mids, hist_near, bin_width)
    axarr[1][0].bar(bin_mids, hist_surf, bin_width)
    axarr[1][1].bar(bin_mids, hist_inte, bin_width)
    axarr[2][0].bar(bin_mids, hist_deep, bin_width)

    for ax1 in axarr:
        for ax2 in ax1:
            ax2.set_xlabel(xtitle)
            ax2.set_ylabel(ytitle)

    axarr[0][0].set_title("Free")
    axarr[0][1].set_title("Near")
    axarr[1][0].set_title("Surface")
    axarr[1][1].set_title("Intermediate")
    axarr[2][0].set_title("Deep")

    axarr[-1][-1].axis('off')

def suq_curve(q, N, A, kc, gamma):
    return (N/A) / (kc*q**4 + gamma*q**2)

def graph_seed_membranemodes(sdlabel,
                             df,
                             nlipids_per_leaflet,
                             ax,
                             mtitle = '',
                             xtitle = '',
                             ytitle = '',
                             color = 'b',
                             xlabel = True):
    r""" Plotting function for total membrane modes
    """
    ax.scatter(df['x_fft'], df['su_fft'], color = 'r', marker = '+', linewidth = 1)
    ax.scatter(df['x_direct'], df['su_direct'], color = 'b', marker = 'o', s = 80, facecolors = 'none')

    x_fft = df['x_fft'].to_numpy()
    su_fft = df['su_fft'].to_numpy()
    qcutoff_mean = np.mean(df['qcutoff_fft'].to_numpy())
    area_mean = np.mean(df['membrane_area'].to_numpy())

    # Try to do a fit
    from scipy.optimize import curve_fit
    # XXX: Hardcoded for now
    popt, pcov = curve_fit(lambda q, kc, gamma: suq_curve(q, nlipids_per_leaflet, area_mean, kc, gamma), x_fft[1:], su_fft[1:], bounds = ([0.0, -np.inf], [np.inf, np.inf]))
    print(f"Simuation fit values:")
    print(f"  kc:    {popt[0]}")
    print(f"  gamma: {popt[1]}")
    ax.plot(x_fft[1:], suq_curve(x_fft[1:], N = nlipids_per_leaflet, A = area_mean, kc = popt[0], gamma = popt[1]), 'r--')

    # Plot the vertical line for qcutoff
    ax.avxline(x = qcutoff_mean, ymin = 0, ymax = 1.0, color = 'm', linestyle = '--')

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_title(mtitle)
    if xlabel: ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

def msd_fit(t, n, D):
    return 2*n*D*t

def graph_seed_simplespheres_msd(sdlabel,
                                 deltatau,
                                 df,
                                 ax,
                                 mtitle = '',
                                 xtitle = '',
                                 ytitle = '',
                                 color = 'b',
                                 xlabel = True):
    r""" Plotting function for simple spehres MSD
    """
    timesteps = df['timestep'].to_numpy()*deltatau
    msd = df['msd_simple'].to_numpy()
    ax.scatter(timesteps, msd, color = 'b', marker = 'o', s = 80, facecolors = 'none')

    # Put a fit on the curve
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(lambda t, D: msd_fit(t, 3, D), timesteps, msd, 1.0)
    ax.plot(timesteps, msd_fit(timesteps, 3, popt[0]), 'r--', linewidth = 1)
    print(f"Diffusion: {popt[0]}")

    ax.set_title(mtitle)
    ax.set_title(mtitle)
    if xlabel: ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    


