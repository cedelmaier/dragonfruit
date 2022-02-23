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

def graph_seed_temperature(sd,
                           ax,
                           color = 'b',
                           xlabel = True):
    r""" Plots the temperature of the seed as measured by kinetic theory
    """
    graph_seed_scatter(sd.label, sd.df["timestep"], sd.df["T"], ax, mtitle = "Temperature", xtitle = "Timestep", ytitle = "Temperature (kT)", color = color)

def graph_seed_pressure(sd,
                        ax,
                        color = 'b',
                        xlabel = True):
    r""" Plots the pressure of the seed as measured by stress tensor
    """
    graph_seed_scatter(sd.label, sd.df["timestep"], sd.df["P"], ax, mtitle = "Pressure", xtitle = "Timestep", ytitle = "Presure (kT/$\sigma^{3}$)", color = color)

def graph_seed_area(sd,
                    ax,
                    color = 'b',
                    xlabel = 'True'):
    r""" Plots the area of the box in the XY dimension
    """
    graph_seed_scatter(sd.label, sd.df["timestep"], sd.df["membrane_area"], ax, mtitle = "Membrane Area", xtitle = "Timestep", ytitle = "Membrane area ($\sigma^{2}$)", color = color)

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

def graph_seed_membranemodes(sd,
                             ax,
                             color = 'b',
                             xlabel = True):
    r""" Plotting function for total membrane modes
    """
    x_fft   = sd.df['x_fft'].to_numpy()
    su_fft  = sd.df['su_fft'].to_numpy()
    x_fft   = x_fft[~np.isnan(x_fft)]
    su_fft  = su_fft[~np.isnan(su_fft)]

    x_direct    = sd.df['x_direct'].to_numpy()
    su_direct   = sd.df['su_direct'].to_numpy()
    x_direct    = x_direct[~np.isnan(x_direct)]
    su_direct   = su_direct[~np.isnan(su_direct)]

    ax.scatter(x_fft[1:], su_fft[1:], color = color, marker = '+', linewidth = 1)
    ax.scatter(x_direct[1:], su_direct[1:], color = color, marker = 'o', s = 80, facecolors = 'none')

    # Figure out where cutoff etc are
    qcutoff_mean = np.nanmean(sd.df['uq_2d_fft_qcutoff'].to_numpy())
    area_mean = np.nanmean(sd.df['membrane_area'].to_numpy())

    # Figure out where the lower cutoff is
    idx = np.where(np.greater(x_fft, qcutoff_mean))
    idx = np.int32(idx[0][0])
    jdx = np.where(np.greater(x_fft, 1.0))
    jdx = np.int32(jdx[0][0])

    # Generate some guesses for the fit
    nlipids_per_leaflet = sd.lipids.nlipids_per_leaflet
    kcguess1 = nlipids_per_leaflet / area_mean / su_fft[idx] / (x_fft[idx]**4)

    # Try to do 2 fits, with and without the surface tension term
    from scipy.optimize import curve_fit
    popt_kc, pcov_kc = curve_fit(lambda q, kc: suq_curve(q, nlipids_per_leaflet, area_mean, kc, 0.0), x_fft[idx:jdx], su_fft[idx:jdx], bounds = ([0.0, np.inf]), p0 = [kcguess1])
    popt_ga, pcov_ga = curve_fit(lambda q, kc, gamma: suq_curve(q, nlipids_per_leaflet, area_mean, kc, gamma), x_fft[idx:jdx], su_fft[idx:jdx], bounds = ([0.0, -np.inf], [np.inf, np.inf]), p0 = [kcguess1, 0.0])

    print(f"Simuation fit values:")
    print(f"  kc(guess):    {kcguess1}")
    print(f"  kc only:      {popt_kc[0]}")
    print(f"  kc (both):    {popt_ga[0]}")
    print(f"  gamma (both): {popt_ga[1]}")
    ax.plot(x_fft[idx:jdx], suq_curve(x_fft[idx:jdx], N = nlipids_per_leaflet, A = area_mean, kc = popt_kc[0], gamma = 0.0), color = color, linestyle = '--')
    ax.plot(x_fft[idx:jdx], suq_curve(x_fft[idx:jdx], N = nlipids_per_leaflet, A = area_mean, kc = popt_ga[0], gamma = popt_ga[1]), color = color, linestyle = ':')

    # Plot the vertical line for qcutoff
    ax.axvline(x = qcutoff_mean, ymin = 0, ymax = 1.0, color = 'k', linestyle = '-')

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_title("Membrane Modes")
    if xlabel: ax.set_xlabel(r'q ($\sigma^{-1}$)')
    ax.set_ylabel(r'$ \langle | u(q) |^{2} \rangle $ ($\sigma^{2}$)')

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
    


