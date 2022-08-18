# Something

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Get a smoothing function
from scipy.signal import savgol_filter

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
from common import moving_average

# General scatter plotting function
def graph_seed_scatter(sdlabel,
                       xdata,
                       ydata,
                       ax,
                       mtitle = '',
                       xtitle = '',
                       ytitle = '',
                       color = 'b',
                       smooth = False,
                       smooth_window = 10,
                       smooth_poly = 3,
                       xlabel = True):
                  
    r""" Generic 1D plotting function for data
    """
    ax.set_title(mtitle)
    ax.set_ylabel(ytitle)
    if xlabel: ax.set_xlabel(xtitle)

    # Check for smoothing options
    ydata_new = ydata.to_numpy().flatten()
    if smooth:
        ydata_new = savgol_filter(ydata, smooth_window, smooth_poly)

    ax.scatter(xdata, ydata_new, color = color, label = sdlabel)

    return ydata_new

# General line plotting function
def graph_seed_plot(sdlabel,
                    xdata,
                    ydata,
                    ax,
                    mtitle = '',
                    xtitle = '',
                    ytitle = '',
                    color = 'b',
                    smooth = False,
                    smooth_window = 10,
                    smooth_poly = 3,
                    xlabel = True):
    r""" Generic 1D plot function
    """
    ax.set_title(mtitle)
    ax.set_ylabel(ytitle)
    if xlabel: ax.set_xlabel(xtitle)

    # Check for smoothing options
    ydata_new = ydata.to_numpy().flatten()
    if smooth:
        ydata_new = savgol_filter(ydata, smooth_window, smooth_poly)

    ax.plot(xdata, ydata_new, color = color, label = sdlabel)

    return ydata_new

# Graph the absolute Z distance
def graph_seed_zdist(sd,
                     ax,
                     color = 'b',
                     xlabel = True):
    r""" Plot the absolute z position of the helix COM
    """
    z = np.abs(sd.master_time_df['helix_z'] - sd.master_time_df['lipid_z'])

    ydata = graph_seed_scatter(sd.label, sd.master_time_df.index/1000.0, z, ax, mtitle = "Absolute Z position", xtitle = "Time (ns)", ytitle = r"Z ($\AA$)", color = color)
    return ydata

# Plot the center of mass distance between lipid and helix, with phosphate groups
def graph_seed_zpos_wheads(sd,
                           ax,
                           color = 'b',
                           xlabel = True):
    r""" Plot the center of mass distance between lipids and helix
    """
    z_protein   = sd.master_time_df['helix_z']
    z_leaf0     = sd.master_time_df['leaflet0_z']
    z_leaf1     = sd.master_time_df['leaflet1_z']
    z_lipid     = sd.master_time_df['lipid_z']
    # Subtract off the position of the lipid COM from everybody else
    z_protein = z_protein - z_lipid
    z_leaf0 = z_leaf0 - z_lipid
    z_leaf1 = z_leaf1 - z_lipid

    #leaf0 = graph_seed_scatter(sd.label, sd.master_time_df.index/1000.0, z_leaf0, ax, mtitle = "", xtitle = "", ytitle = r"", color = 'k', smooth = True, smooth_window = 11, smooth_poly = 3)
    #leaf1 = graph_seed_scatter(sd.label, sd.master_time_df.index/1000.0, z_leaf1, ax, mtitle = "", xtitle = "", ytitle = r"", color = 'k', smooth = True, smooth_window = 11, smooth_poly = 3)
    #ydata = graph_seed_scatter(sd.label, sd.master_time_df.index/1000.0, z_protein, ax, mtitle = "Z-distance", xtitle = "Time (ns)", ytitle = r"Z ($\AA$)", color = color, smooth = True, smooth_window = 11, smooth_poly = 3)
    leaf0 = graph_seed_plot(sd.label, sd.master_time_df.index/1000.0, z_leaf0, ax, mtitle = "", xtitle = "", ytitle = r"", color = 'k', smooth = True, smooth_window = 11, smooth_poly = 3)
    leaf1 = graph_seed_plot(sd.label, sd.master_time_df.index/1000.0, z_leaf1, ax, mtitle = "", xtitle = "", ytitle = r"", color = 'k', smooth = True, smooth_window = 11, smooth_poly = 3)
    ydata = graph_seed_plot(sd.label, sd.master_time_df.index/1000.0, z_protein, ax, mtitle = "Z-distance", xtitle = "Time (ns)", ytitle = r"Z ($\AA$)", color = color, smooth = True, smooth_window = 11, smooth_poly = 3)

    # Set the appropriate limits
    ylow    = -30.0
    yhi     = 30.0
    ax.set_ylim([ylow, yhi])

    return [ylow, yhi, [leaf0, leaf1, ydata]]

# Graph the fracitonal helicity
def graph_seed_helicity(sd,
                        ax,
                        color = 'b',
                        xlabel = True):
    r""" Plot the fractional helicity of the helix
    """
    helicity = sd.master_time_df['helicity']

    #graph_seed_scatter(sd.label, sd.master_time_df.index/1000.0, helicity, ax, mtitle = "Helicity", xtitle = "Time (ns)", ytitle = r"Helicity (AU)", color = color)
    ydata = graph_seed_plot(sd.label, sd.master_time_df.index/1000.0, helicity, ax, mtitle = "Helicity", xtitle = "Time (ns)", ytitle = r"Helicity (AU)", color = color, smooth = True, smooth_window = 11, smooth_poly = 3)

    # Set the appropriate limits
    ylow    = 0.0
    yhi     = 1.05
    ax.set_ylim([ylow, yhi])

    return [ylow, yhi, ydata]

# Graph the global tilt
def graph_seed_globaltilt(sd,
                          ax,
                          color = 'b',
                          xlabel = True):
    r""" Plot the global tilt (as determined by HELANAL)
    """
    global_tilt = sd.master_time_df['global_tilt']

    #ydata = graph_seed_scatter(sd.label, sd.master_time_df.index/1000.0, global_tilt, ax, mtitle = "Global Tilt", xtitle = "Time (ns)", ytitle = r"Global Tilt (deg)", color = color, smooth = True, smooth_window = 11, smooth_poly = 3)
    ydata = graph_seed_plot(sd.label, sd.master_time_df.index/1000.0, global_tilt, ax, mtitle = "Global Tilt", xtitle = "Time (ns)", ytitle = r"Global Tilt (deg)", color = color, smooth = True, smooth_window = 11, smooth_poly = 3)

    # Set the appropriate limits
    ylow    = 0.0
    yhi     = 180.0
    ax.set_ylim([ylow, yhi])

    return [ylow, yhi, ydata]
