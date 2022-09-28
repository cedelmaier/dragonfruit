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
                       xlabel = True,
                       alpha = 1.0):
                  
    r""" Generic 1D plotting function for data
    """
    ax.set_title(mtitle)
    ax.set_ylabel(ytitle)
    if xlabel: ax.set_xlabel(xtitle)

    # Check for smoothing options
    ydata_new = ydata.to_numpy().flatten()
    if smooth:
        ydata_new = savgol_filter(ydata, smooth_window, smooth_poly)

    ax.scatter(xdata, ydata_new, color = color, label = sdlabel, alpha = alpha)

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
                    xlabel = True,
                    alpha = 1.0,
                    linestyle = 'solid'):
    r""" Generic 1D plot function
    """
    ax.set_title(mtitle)
    ax.set_ylabel(ytitle)
    if xlabel: ax.set_xlabel(xtitle)

    # Check for smoothing options
    ydata_new = ydata.to_numpy().flatten()
    if smooth:
        ydata_new = savgol_filter(ydata, smooth_window, smooth_poly)

    ax.plot(xdata, ydata_new, color = color, label = sdlabel, alpha = alpha, linestyle = linestyle)

    return ydata_new

# Graph the absolute Z distance
def graph_seed_zdist(sd,
                     ax,
                     color = 'b',
                     xlabel = True,
                     dolegend = False):
    r""" Plot the absolute z position of the helix COM
    """
    z = np.abs(sd.master_time_df['helix_z'] - sd.master_time_df['lipid_z'])

    ydata = graph_seed_scatter(sd.label, sd.master_time_df.index/1000.0, z, ax, mtitle = "Absolute Z position", xtitle = "Time (ns)", ytitle = r"Z ($\AA$)", color = color)
    return ydata

# Plot the center of mass distance between lipid and helix, with phosphate groups
def graph_seed_zpos_wheads(sd,
                           ax,
                           color = 'b',
                           xlabel = True,
                           dolegend = False):
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
                        xlabel = True,
                        dolegend = False):
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
                          xlabel = True,
                          dolegend = False):
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

# Graph the pdipole tilt
def graph_seed_pdipoletilt(sd,
                            ax,
                            color = 'b',
                            xlabel = True,
                            dolegend = False):
    r""" Plot the pdipole tilt (Z-axis as reference)
    """
    p_dipole = sd.master_time_df[['p_dipole_x', 'p_dipole_y', 'p_dipole_z']].to_numpy()
    z_axis = np.array([0.0, 0.0, 1.0])
    ntimes = len(sd.master_time_df.index)
    angles = np.zeros(ntimes)
    for idx in range(ntimes):
        angle = np.arccos(np.dot(p_dipole[idx,:], z_axis) / (np.linalg.norm(p_dipole[idx,:]) * np.linalg.norm(z_axis)))
        angle = angle * 180.0/np.pi
        angles[idx] = angle

    # Have to do some shennanigans to get the correct form for our plotting
    angles_series = pd.Series(angles)
    ydata = graph_seed_plot(sd.label, sd.master_time_df.index/1000.0, angles_series, ax, mtitle = "Helix Dipole Tilt", xtitle = "Time (ns)", ytitle = r"Helix Dipole Tilt (deg)", color = color, smooth = True, smooth_window = 11, smooth_poly = 3)

    # Set the appropriate limits
    ylow    = 0.0
    yhi     = 180.0
    ax.set_ylim([ylow, yhi])

    return [ylow, yhi, ydata]

# Graph the pdipole tilt
def graph_seed_helixpdipoleangle(sd,
                            ax,
                            color = 'b',
                            xlabel = True,
                            dolegend = False):
    r""" Plot the angle between the helix and p_dipole
    """
    helix_axis = sd.master_time_df[['helix_global_axis_x', 'helix_global_axis_y', 'helix_global_axis_z']].to_numpy()
    p_dipole = sd.master_time_df[['p_dipole_x', 'p_dipole_y', 'p_dipole_z']].to_numpy()
    ntimes = len(sd.master_time_df.index)
    angles = np.zeros(ntimes)
    for idx in range(ntimes):
        angle = np.arccos(np.dot(p_dipole[idx,:], helix_axis[idx,:]) / (np.linalg.norm(p_dipole[idx,:]) * np.linalg.norm(helix_axis[idx,:])))
        angle = angle * 180.0/np.pi
        angles[idx] = angle

    # Have to do some shennanigans to get the correct form for our plotting
    angles_series = pd.Series(angles)
    ydata = graph_seed_plot(sd.label, sd.master_time_df.index/1000.0, angles_series, ax, mtitle = "Helix-Dipole Angle", xtitle = "Time (ns)", ytitle = r"Helix-Dipole Angle (deg)", color = color, smooth = True, smooth_window = 11, smooth_poly = 3)

    # Set the appropriate limits
    ylow    = 0.0
    yhi     = 180.0
    ax.set_ylim([ylow, yhi])

    return [ylow, yhi, ydata]

# Graph the pdipole moment
def graph_seed_pdipolemoment(sd,
                            ax,
                            color = 'b',
                            xlabel = True,
                            dolegend = False):
    r""" Plot the electric dipole moment of helix
    """
    p_dipole = sd.master_time_df[['p_dipole_x', 'p_dipole_y', 'p_dipole_z']].to_numpy()
    ntimes = len(sd.master_time_df.index)
    magnitudes = np.zeros(ntimes)
    for idx in range(ntimes):
        mag = np.linalg.norm(p_dipole[idx,:])
        magnitudes[idx] = mag

    # Have to do some shennanigans to get the correct form for our plotting
    mag_series = pd.Series(magnitudes)
    ydata = graph_seed_plot(sd.label, sd.master_time_df.index/1000.0, mag_series, ax, mtitle = "Helix Electric Dipole Magnitude", xtitle = "Time (ns)", ytitle = r"Helix Electric Dipole Magnitude", color = color, smooth = True, smooth_window = 11, smooth_poly = 3)

    # Set the appropriate limits
    ylow    = 0.0
    yhi     = 40.0
    ax.set_ylim([ylow, yhi])

    return [ylow, yhi, ydata]

# Graph the Z-component of the force on the helix
def graph_seed_zforce(sd,
                      ax,
                      color = 'b',
                      xlabel = True,
                      dolegend = True):
    r""" Plot the Z-component of the forces for all, electrostatic, and other contributions
    """
    # Keep the force types in this order
    force_types = ["all", "q", "noq"]
    #alphas = [1.0, 2.0/3.0, 1.0/3.0]
    alphas = [1.0, 1.0, 1.0]
    linestyles = ["solid", "dotted", "dashed"]
    ntimes = len(sd.master_forces_df.index)
    ydata_common = []
    for (force_type,alpha,linestyle) in zip(force_types, alphas, linestyles):
        target_name = "force_com_" + force_type + "_z"
        zf = sd.master_forces_df[[target_name]].to_numpy()
        zforce = np.zeros(ntimes)
        for idx in range(ntimes):
            zforce[idx] = zf[idx]
        z_series = pd.Series(zforce)
        ydata = graph_seed_plot(force_type, sd.master_forces_df.index/1000.0, z_series, ax, mtitle = "Helix Force Z-component", xtitle = "Time (ns)", ytitle = r"Force (kJ mol$^{-1}$ nm$^{-1}$)", color = color, alpha = alpha, linestyle = linestyle, smooth = True, smooth_window = 11, smooth_poly = 3)
        ydata_common.append(ydata)
    if dolegend: ax.legend()

    # Set the appropriate limits
    ylow    = -8000.0
    yhi     = 8000.0
    ax.set_ylim([ylow, yhi])

    return [ylow, yhi, ydata_common]

# Graph the torque perpendicular to the uhat and z-axis
def graph_seed_perptorque(sd,
                          ax,
                          color = 'b',
                          xlabel = True,
                          dolegend = True):
    r""" Plot the torque component perpendicular to the helix axis (uhat) and z-axis (zhat)
    """
    # Break down forces into components
    force_types = ["all", "q", "noq"]
    #alphas = [1.0, 2.0/3.0, 1.0/3.0]
    alphas = [1.0, 1.0, 1.0]
    linestyles = ["solid", "dotted", "dashed"]

    # Get the two orientation we need, however, relaize that they have different timepoints
    zhat = np.array([0.0, 0.0, 1.0])
    uhat_df = sd.master_time_df[['helix_global_axis_x', 'helix_global_axis_y', 'helix_global_axis_z']]
    ntimes_forces = len(sd.master_forces_df.index)
    ydata_common = []
    for (force_type,alpha,linestyle) in zip(force_types, alphas, linestyles):
        target_name = "torque_com_" + force_type + "_*"
        torque_df_filtered = sd.master_forces_df.filter(regex = target_name)
        uhat_df_filtered = uhat_df[uhat_df.index.isin(torque_df_filtered.index)]

        # Now that we have our dataframes, we can construct out variables
        uhat_arr    = uhat_df_filtered.to_numpy()
        torque_arr  = torque_df_filtered.to_numpy()

        torque_dot_nhat = np.zeros(ntimes_forces)
        for idx in range(ntimes_forces):
            # Construct the normalized uhat
            uhat = uhat_arr[idx,:] / np.linalg.norm(uhat_arr[idx,:])
            nhat = np.cross(uhat, zhat)
            nhat /= np.linalg.norm(nhat)

            torque_dot_nhat[idx] = np.dot(torque_arr[idx,:], nhat)

        # Convert to a series and plot
        torque_series = pd.Series(torque_dot_nhat)
        ydata = graph_seed_plot(force_type, sd.master_forces_df.index/1000.0, torque_series, ax, mtitle = "Helix Perpendicular Torque", xtitle = "Time (ns)", ytitle = r"Torque (kJ mol$^{-1}$)", color = color, alpha = alpha, linestyle = linestyle, smooth = True, smooth_window = 11, smooth_poly = 3)
        ydata_common.append(ydata)
    if dolegend: ax.legend()

    # Set the appropriate limits
    ylow    = -800000.0
    yhi     = 800000.0
    ax.set_ylim([ylow, yhi])

    return [ylow, yhi, ydata_common]
