#!/usr/bin/env python3

# XXX: Put a license here

""" Simple script to combine the gromacs analyses listed, keep updating and clean up in future """

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

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src', 'lib'))
from stylelib.common_styles import *

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data')

# Wrap up stuff into simulation packages
simulation_packages = {
                       'polarmonomer_aglipid_rotxN_50mMKCl':  [[
                                                              'unbiased_zdepth00_rotx0_helix_50mMKCl',
                                                              'unbiased_zdepth00_rotx90_helix_50mMKCl',
                                                              'unbiased_zdepth00_rotx180_helix_50mMKCl',
                                                              'unbiased_zdepth00_rotx270_helix_50mMKCl',
                                                              ],
                                                              [
                                                              r'AGMAGL 0 deg 50 mM KCl',
                                                              r'AGMAGL 90 deg 50 mM KCl',
                                                              r'AGMAGL 180 deg 50 mM KCl',
                                                              r'AGMAGL 270 deg 50 mM KCl',
                                                              ]],
                       'polarmonomer_aglipid_rotxNtrans_50mMKCl':  [[
                                                              'unbiased_zdepth00_rotx0_helix_50mMKCl',
                                                              'unbiased_zdepth00_rotx90_helix_50mMKCl',
                                                              'unbiased_zdepth00_rotx180_helix_50mMKCl',
                                                              'unbiased_zdepth00_rotx270_helix_50mMKCl',
                                                              'unbiased_zdepth00_trans_helix_50mMKCl',
                                                              ],
                                                              [
                                                              r'AGMAGL 0 deg 50 mM KCl',
                                                              r'AGMAGL 90 deg 50 mM KCl',
                                                              r'AGMAGL 180 deg 50 mM KCl',
                                                              r'AGMAGL 270 deg 50 mM KCl',
                                                              r'AGMAGL Trans 50 mM KCl',
                                                              ]],
                       'polarmonomer_aglipid_rotxN_150mMKCl':  [[
                                                              'unbiased_150mMKCl',
                                                              'agmonomer_11x11_zdepth00_rotx90_150KCl',
                                                              'agmonomer_aglipid_11x11_zdepth00_rotx180_150mMKCl',
                                                              'agmonomer_aglipid_11x11_zdepth00_rotx270_150mMKCl',
                                                              ],
                                                              [
                                                              r'AGMAGL 0 deg 150 mM KCl',
                                                              r'AGMAGL 90 deg 150 mM KCl',
                                                              r'AGMAGL 180 deg 150 mM KCl',
                                                              r'AGMAGL 270 deg 150 mM KCl',
                                                              ]],
                       'neutralmonomer_aglipid_rotxN_50mMKCl':  [[
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx0_50mMKCl',
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx90_50mMKCl',
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx175_50mMKCl',
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx270_50mMKCl',
                                                              ],
                                                              [
                                                              r'RFMAGL 0 deg 50 mM KCl',
                                                              r'RFMAGL 90 deg 50 mM KCl',
                                                              r'RFMAGL 175 deg 50 mM KCl',
                                                              r'RFMAGL 270 deg 50 mM KCl',
                                                              ]],
                       'neutralmonomer_aglipid_rotxNall_50mMKCl':  [[
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx0_50mMKCl',
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx90_50mMKCl',
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx175_50mMKCl',
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx265_50mMKCl',
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx270_50mMKCl',
                                                              ],
                                                              [
                                                              r'RFMAGL 0 deg 50 mM KCl',
                                                              r'RFMAGL 90 deg 50 mM KCl',
                                                              r'RFMAGL 175 deg 50 mM KCl',
                                                              r'RFMAGL 265 deg 50 mM KCl',
                                                              r'RFMAGL 270 deg 50 mM KCl',
                                                              ]],
                        'neutralmonomer_aglipid_rotxN_150mMKCl': [[
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx0_150KCl',
                                                              ],
                                                              [
                                                              r'RFMAGL 0 deg 150 mM KCl',
                                                              ]],
                                                              }

# Get the simulations to run
simulations_to_run = [
                      'polarmonomer_aglipid_rotxN_50mMKCl',
                      #'polarmonomer_aglipid_rotxNtrans_50mMKCl',
                      'polarmonomer_aglipid_rotxN_150mMKCl',
                      'neutralmonomer_aglipid_rotxN_50mMKCl',
                      #'neutralmonomer_aglipid_rotxNall_50mMKCl',
                      'neutralmonomer_aglipid_rotxN_150mMKCl',
                      ]
simulation_text_data = {}
simulation_text_data['polarmonomer_aglipid_rotxN_50mMKCl'] = 'Polar 50 mM KCl'
simulation_text_data['polarmonomer_aglipid_rotxNtrans_50mMKCl'] = 'Polar 50 mM KCl (with trans)'
simulation_text_data['polarmonomer_aglipid_rotxN_150mMKCl'] = 'Polar 150 mM KCl'
simulation_text_data['neutralmonomer_aglipid_rotxN_50mMKCl'] = 'Neutral 50 mM KCl'
simulation_text_data['neutralmonomer_aglipid_rotxNall_50mMKCl'] = 'Neutral 50 mM KCl (all)'
simulation_text_data['neutralmonomer_aglipid_rotxN_150mMKCl'] = 'Neutral 150 mM KCl'

# Some final resting places for data, last 25 ns of the simulations
measurement_names = [
                     'zpos',
                     'helix',
                     'tilt',
                     'pdip',
                     'hp',
                    ]
measurement_text = {}
measurement_text['zpos'] = (r"Z ($\AA$)", -30.0, 30.0)
measurement_text['helix'] = (r"Helicity (AU)", 0.0, 1.05)
measurement_text['tilt'] = (r"Helix tilt (deg)", 0.0, 180.0)
measurement_text['pdip'] = (r"Helix dipole tilt (deg)", 0.0, 180.0)
measurement_text['hp'] = (r"Helix-dipole angle (deg)", 0.0, 180.0)

final_data_dict = {}
final_data      = {}
for measurement_name in measurement_names:
    final_data_dict[measurement_name] = {}
    final_data[measurement_name] = {}

# Common parameter choices for end times, things like that
end_times = 175000.0

for simtorun in simulations_to_run:
    print(f"Simulation block {simtorun}")

    # Set up the plots beforehand
    #plt.style.use(septin_poster_stl)
    plt.style.use(septin_runs_stl)
    fig_zdist, ax_zdist = plt.subplots(1, 1, figsize = (15, 10))
    fig_zpos , ax_zpos  = plt.subplots(1, 1, figsize = (15, 10))
    fig_twist, ax_twist = plt.subplots(1, 1, figsize = (15, 10))
    fig_resid, ax_resid = plt.subplots(1, 1, figsize = (15, 10))
    fig_helix, ax_helix = plt.subplots(1, 1, figsize = (15, 10))
    fig_tilt , ax_tilt  = plt.subplots(1, 1, figsize = (15, 10))    # Helix tilt
    fig_pdip , ax_pdip  = plt.subplots(1, 1, figsize = (15, 10))    # Dipole tilt
    fig_hp   , ax_hp    = plt.subplots(1, 1, figsize = (15, 10))    # Helix-dipole angle
    
    # Set up a stride and only graph every stride points
    stride = 10
    
    # For each simulation, loop through and get the relevant data and plot it
    Nsims = 4
    cidx = 0
    current_simulation_parameters = simulation_packages[simtorun]

    outname = simtorun
    simulation_names = current_simulation_parameters[0][:]
    simulation_legends = current_simulation_parameters[1][:]

    # Safety check if we ever want to reuse data!
    for k,v in final_data_dict.items():
        if simtorun not in v:
            v[simtorun] = np.array([], dtype=np.float64)

    for isim,sname in enumerate(simulation_names):
        for iseed in np.arange(1, Nsims+1):
            print(f"Sim: {isim} = {sname}, N{iseed}")
            filepath = os.path.join(main_path, sname, f"N{iseed}")
    
            hd5_filename = os.path.join(filepath, f"N{iseed}.h5")
            master_time_df = pd.read_hdf(hd5_filename)
    
            times = master_time_df.index
            flat_times = times.to_numpy().flatten()
    
            labelname = simulation_legends[isim] + f" N{iseed}"
    
            # Compute the zdist plot
            zdist = np.abs(master_time_df['helix_z'] - master_time_df['lipid_z'])
            ax_zdist.plot(times, zdist, label = labelname)
    
            # Compute the z position
            z_pos           = master_time_df[['helix_z']].to_numpy().flatten()
            leaflet0_pos    = master_time_df[['leaflet0_z']].to_numpy().flatten()
            leaflet1_pos    = master_time_df[['leaflet1_z']].to_numpy().flatten()
            lipid_pos       = master_time_df[['lipid_z']].to_numpy().flatten()
            z_pbc           = master_time_df[['unit_cell_z']].to_numpy().flatten()
            # Subtract off the lipid COM position
            z_pos = z_pos - lipid_pos
            leaflet0_pos = leaflet0_pos - lipid_pos
            leaflet1_pos = leaflet1_pos - lipid_pos
            z_pbc = z_pbc/2 - lipid_pos
            # Correct the position if under the lower leaflet
            for idx in range(len(z_pos)):
                if z_pos[idx] < leaflet1_pos[idx]:
                    z_pos[idx] = z_pbc[idx] - z_pos[idx]
    
            # Plot the location of everything at some stride distance
            ax_zpos.plot(times[::stride], z_pos[::stride], label = labelname)
            ax_zpos.plot(times[::stride], leaflet0_pos[::stride], color = 'k')
            ax_zpos.plot(times[::stride], leaflet1_pos[::stride], color = 'k')
    
            # Compute the helicity
            helicity = master_time_df[['helicity']].to_numpy().flatten()
            # Smooth this data, as it is gross right now
            ax_helix.plot(times[::stride], helicity[::stride], label = labelname)
    
            # Get the tilt
            global_tilt = master_time_df[['global_tilt']].to_numpy().flatten()
            ax_tilt.plot(times[::stride], global_tilt[::stride], label = labelname)
    
            # Get and calculate the dipole tilt and helix-dipole angle
            ntimes = len(master_time_df.index)
            pdipole_tilt_angles     = np.zeros(ntimes)
            helix_pdipole_angles    = np.zeros(ntimes)
            pdipole_arr             = master_time_df[['p_dipole_x', 'p_dipole_y', 'p_dipole_z']].to_numpy()
            helix_arr               = master_time_df[['helix_global_axis_x', 'helix_global_axis_y', 'helix_global_axis_z']].to_numpy()
            z_axis                  = np.array([0.0, 0.0, 1.0])
            for itime in range(ntimes):
                pdipole_angle = np.arccos(np.dot(pdipole_arr[itime,:], z_axis) / (np.linalg.norm(pdipole_arr[itime,:]) * np.linalg.norm(z_axis)))
                helixpd_angle = np.arccos(np.dot(pdipole_arr[itime,:], helix_arr[itime,:]) / (np.linalg.norm(pdipole_arr[itime,:]) * np.linalg.norm(helix_arr[itime,:])))
                pdipole_tilt_angles[itime] = 180.0 / np.pi * pdipole_angle
                helix_pdipole_angles[itime] = 180.0 / np.pi * helixpd_angle
            ax_pdip.plot(times[::stride], pdipole_tilt_angles[::stride], label = labelname)
            ax_hp.plot(times[::stride], helix_pdipole_angles[::stride], label = labelname)

            # Get the final values of some things and record them to process later
            # Do some fancy boolean indexing on the times
            final_indices = np.where(flat_times < end_times)
            final_data_dict['zpos'][simtorun]   = np.concatenate((final_data_dict['zpos'][simtorun],    z_pos[final_indices]), axis = 0)
            final_data_dict['helix'][simtorun]  = np.concatenate((final_data_dict['helix'][simtorun],   helicity[final_indices]), axis = 0)
            final_data_dict['tilt'][simtorun]   = np.concatenate((final_data_dict['tilt'][simtorun],    global_tilt[final_indices]), axis = 0)
            final_data_dict['pdip'][simtorun]   = np.concatenate((final_data_dict['pdip'][simtorun],    pdipole_tilt_angles[final_indices]), axis = 0)
            final_data_dict['hp'][simtorun]     = np.concatenate((final_data_dict['hp'][simtorun],      helix_pdipole_angles[final_indices]), axis = 0)

    
        cidx += 1

    ax_zdist.set_xlabel('Time (ps)')
    ax_zdist.set_ylabel(r'Z distance ($\AA$)')
    #ax_zdist.legend(loc = 'lower right')
    fig_zdist.tight_layout()
    fig_zdist.savefig('gromacs_zdist_' + outname + '.pdf', dpi = fig_zdist.dpi)
    plt.close()
    
    ax_zpos.set_xlabel('Time (ps)')
    ax_zpos.set_ylabel(r'z ($\AA$)')
    #ax_zpos.legend(loc = 'lower right')
    fig_zpos.tight_layout()
    fig_zpos.savefig('gromacs_zpos_' + outname + '.pdf', dpi = fig_zpos.dpi)
    plt.close()
    
    ax_helix.set_xlabel('Time (ps)')
    ax_helix.set_ylabel('Helical nature (AU)')
    #ax_helix.legend(loc = 'lower right')
    ax_helix.set_ylim([0.0, 1.1])
    fig_helix.tight_layout()
    fig_helix.savefig('gromacs_helicity_' + outname + '.pdf', dpi = fig_helix.dpi)
    plt.close()
    
    ax_tilt.set_xlabel('Time (ps)')
    ax_tilt.set_ylabel('Global tilt (deg)')
    #ax_tilt.legend(loc = 'lower right')
    ax_tilt.set_ylim([0.0, 180.0])
    fig_tilt.tight_layout()
    fig_tilt.savefig('gromacs_tilt_' + outname + '.pdf', dpi = fig_tilt.dpi)
    plt.close()
    
    ax_pdip.set_xlabel('Time (ps)')
    ax_pdip.set_ylabel('Helix dipole tilt (deg)')
    #ax_tilt.legend(loc = 'lower right')
    ax_pdip.set_ylim([0.0, 180.0])
    fig_pdip.tight_layout()
    fig_pdip.savefig('gromacs_pdipoletilt_' + outname + '.pdf', dpi = fig_pdip.dpi)
    plt.close()
    
    ax_hp.set_xlabel('Time (ps)')
    ax_hp.set_ylabel('Helix dipole - Helix axis angle (deg)')
    #ax_tilt.legend(loc = 'lower right')
    ax_hp.set_ylim([0.0, 180.0])
    fig_hp.tight_layout()
    fig_hp.savefig('gromacs_helixpdipoleangle_' + outname + '.pdf', dpi = fig_hp.dpi)
    plt.close()

    # Calculate some final mean/std variables
    for measurement_name in measurement_names:
        final_data[measurement_name][simtorun] = (np.mean(final_data_dict[measurement_name][simtorun]),
                                                  np.std(final_data_dict[measurement_name][simtorun], ddof=1),
                                                  np.std(final_data_dict[measurement_name][simtorun], ddof=1))

#  Plot things like the final measurements
for measurement_name in measurement_names:
    print(f"Measurement name: {measurement_name}, results: {final_data[measurement_name]}")
    m_fig, m_ax = plt.subplots(1, 1, figsize = (15, 10))
    xvals = []
    yvals = []
    yerrs = []
    for simtorun in simulations_to_run:
        xvals.append(simulation_text_data[simtorun])
        yvals.append(final_data[measurement_name][simtorun][0])
        yerrs.append(final_data[measurement_name][simtorun][1])
    print(f"  {yvals} +/- {yerrs}")
    m_ax.scatter(xvals, yvals, zorder = 100,
                 s = 100, marker = 's', color = 'k', label = None)
    m_ax.errorbar(xvals, yvals, yerr = yerrs,
                  ecolor = 'k', elinewidth = 2, capsize = 7, capthick = 1, zorder = 0,
                  fmt = 'none')
    m_ax.set_ylabel(measurement_text[measurement_name][0])
    m_ax.set_ylim([measurement_text[measurement_name][1], measurement_text[measurement_name][2]])
    m_fig.tight_layout()
    m_fig.savefig('final_' + measurement_name + '.pdf', dpi = m_fig.dpi)
    plt.close()
