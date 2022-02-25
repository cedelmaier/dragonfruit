# Something

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from common import moving_average

def suq_curve(q, N, A, kc, gamma):
    return (N/A) / (kc*q**4 + gamma*q**2)

def graph_sim_membranemodes(sim,
                            ax,
                            color = 'b',
                            **kwargs):
    r""" Plotting function for membrane modes average over several simulations
    """
    # Prepare some variables from the first seed
    min_size_fft    = 0
    #min_size_direct = 0
    qcutoff_mean = []
    area_mean = []
    for sd in sim.seeds:
        fft_length = sd.df['x_fft'].to_numpy()
        fft_length = fft_length[~np.isnan(fft_length)]
        fft_length = len(fft_length)
        if fft_length > min_size_fft: min_size_fft = fft_length

        #direct_length = sd.df['x_direct'].to_numpy()
        #direct_length = direct_length[~np.isnan(direct_length)]
        #direct_length = len(direct_length)
        #if direct_length > min_size_direct: min_size_direct = direct_length

        qcutoff_mean.append(np.nanmean(sd.df['uq_2d_fft_qcutoff'].to_numpy()))
        area_mean.append(np.nanmean(sd.df['membrane_area'].to_numpy()))

    qcutoff_mean = np.mean(qcutoff_mean)
    area_mean = np.mean(area_mean)

    x_fft_arr = None
    su_fft_arr = None
    #x_direct_arr = None
    #su_direct_arr = None

    for sd in sim.seeds:
        # Get the data for each seed and add it to an average
        nlipids_per_leaflet = sd.lipids.nlipids_per_leaflet

        x_fft   = sd.df['x_fft'].to_numpy()
        su_fft  = sd.df['su_fft'].to_numpy()
        x_fft   = x_fft[~np.isnan(x_fft)]
        su_fft  = su_fft[~np.isnan(su_fft)]

        #x_direct    = sd.df['x_direct'].to_numpy()
        #su_direct   = sd.df['su_direct'].to_numpy()
        #x_direct    = x_direct[~np.isnan(x_direct)]
        #su_direct   = su_direct[~np.isnan(su_direct)]

        # Cheap trick to assign the right size of array the first time around
        if type(x_fft_arr) != type(x_fft):
            x_fft_arr       = x_fft
            su_fft_arr      = su_fft
            #x_direct_arr    = x_direct
            #su_direct_arr   = su_direct
        else:
            x_fft_arr       = np.column_stack((x_fft_arr, x_fft))
            su_fft_arr      = np.column_stack((su_fft_arr, su_fft))
            #x_direct_arr    = np.column_stack((x_direct_arr, x_direct))
            #su_direct_arr   = np.column_stack((su_direct_arr, su_direct))

    # Create the average/mean/std variables
    if type(x_fft_arr) == None:
        print(f"Didn't have any FFT components")
    else:
        # X coordinates should not change, but do the average anyway
        x_fft_mean = np.mean(x_fft_arr, axis = 1)
        su_fft_mean = np.mean(su_fft_arr, axis = 1)
        su_fft_std = np.std(su_fft_arr, axis = 1)

        #x_direct_mean = np.mean(x_direct_arr, axis = 1)
        #su_direct_mean = np.mean(su_direct_arr, axis = 1)
        #su_direct_std = np.std(su_direct_arr, axis = 1)

    ax.errorbar(x_fft_mean, su_fft_mean, yerr = su_fft_std, label = 'FFT', color = color, marker = '+', linestyle = 'none')
    #ax.errorbar(x_direct_mean, su_direct_mean, yerr = su_direct_std, label = 'Direct', color = color, marker = 'o', linestyle = 'none')

    # Try to fit the FFT data with the appropriate curves

    idx = np.where(np.greater(x_fft_mean, qcutoff_mean))
    idx = np.int32(idx[0][0])
    jdx = np.where(np.greater(x_fft_mean, 1.0))
    jdx = np.int32(jdx[0][0])
    # Generate some guesses for the fit
    kcguess1 = nlipids_per_leaflet / area_mean / su_fft_mean[idx] / (x_fft_mean[idx]**4)

    # Try to do a fit on this
    from scipy.optimize import curve_fit
    popt_kc, pcov_kc = curve_fit(lambda q, kc: suq_curve(q, nlipids_per_leaflet, area_mean, kc, 0.0), x_fft_mean[idx:jdx], su_fft_mean[idx:jdx], bounds = ([0.0, np.inf]), p0 = [kcguess1])
    popt_ga, pcov_ga = curve_fit(lambda q, kc, gamma: suq_curve(q, nlipids_per_leaflet, area_mean, kc, gamma), x_fft_mean[idx:jdx], su_fft_mean[idx:jdx], bounds = ([0.0, -np.inf], [np.inf, np.inf]), p0 = [kcguess1, 0.0])

    print(f"Simuation fit values:")
    print(f"  kc(guess):    {kcguess1}")
    print(f"  kc only:      {popt_kc[0]}")
    print(f"  kc (both):    {popt_ga[0]}")
    print(f"  gamma (both): {popt_ga[1]}")
    ax.plot(x_fft_mean[idx:jdx], suq_curve(x_fft[idx:jdx], N = nlipids_per_leaflet, A = area_mean, kc = popt_kc[0], gamma = 0.0), color = color, linestyle = '--')
    ax.plot(x_fft[idx:jdx], suq_curve(x_fft_mean[idx:jdx], N = nlipids_per_leaflet, A = area_mean, kc = popt_ga[0], gamma = popt_ga[1]), color = color, linestyle = ':')

    # Set the correct log scale stuff
    ax.set_yscale('log')
    ax.set_xscale('log')

    if 'xlabel' in kwargs and kwargs['xlabel'] == True: ax.set_xlabel(r'q ($\sigma^{-1}$)')


