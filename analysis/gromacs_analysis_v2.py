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

# Set the style for the plots
plt.style.use(septin_poster_stl)

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
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx90_150mMKCl',
                                                              'rfmonomer_aglipid_11x11_zdepth00_rotx270_150mMKCl',
                                                              ],
                                                              [
                                                              r'RFMAGL 0 deg 150 mM KCl',
                                                              r'RFMAGL 90 deg 150 mM KCl',
                                                              r'RFMAGL 270 deg 150 mM KCl',
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
                     'pmom',
                     'hp',
                    ]
measurement_text = {}
measurement_text['zpos'] = (r"Z ($\AA$)", -30.0, 30.0)
measurement_text['helix'] = (r"Helicity (AU)", 0.0, 1.05)
measurement_text['tilt'] = (r"Helix tilt (deg)", 0.0, 180.0)
measurement_text['pdip'] = (r"Helix dipole tilt (deg)", 0.0, 180.0)
measurement_text['pmom'] = (r"Helix electric dipole moment (?)", 0.0, 50.0)
measurement_text['hp'] = (r"Helix-dipole angle (deg)", 0.0, 180.0)

final_data_dict = {}
final_data      = {}
for measurement_name in measurement_names:
    final_data_dict[measurement_name] = {}
    final_data[measurement_name] = {}

# Categorize into tilted/non-tilted
tilt_criteria = 50.0   # Should be 0 to 40, 140 to 180
alpha_plot = 0.25
fig_zpos_tilt, ax_zpos_tilt     = plt.subplots(1, 1, figsize = (15, 10))
fig_helix_tilt, ax_helix_tilt   = plt.subplots(1, 1, figsize = (15, 10))
fig_tilt_tilt, ax_tilt_tilt     = plt.subplots(1, 1, figsize = (15, 10))
fig_pdip_tilt, ax_pdip_tilt     = plt.subplots(1, 1, figsize = (15, 10))
fig_hp_tilt, ax_hp_tilt         = plt.subplots(1, 1, figsize = (15, 10))
fig_pmom_tilt, ax_pmom_tilt     = plt.subplots(1, 1, figsize = (15, 10))

# Common parameter choices for end times, things like that
end_times = 175000.0

for simtorun in simulations_to_run:
    print(f"Simulation block {simtorun}")

    # Set up the plots beforehand
    fig_zdist, ax_zdist = plt.subplots(1, 1, figsize = (15, 10))
    fig_zpos , ax_zpos  = plt.subplots(1, 1, figsize = (15, 10))
    fig_twist, ax_twist = plt.subplots(1, 1, figsize = (15, 10))
    fig_resid, ax_resid = plt.subplots(1, 1, figsize = (15, 10))
    fig_helix, ax_helix = plt.subplots(1, 1, figsize = (15, 10))
    fig_tilt , ax_tilt  = plt.subplots(1, 1, figsize = (15, 10))    # Helix tilt
    fig_pdip , ax_pdip  = plt.subplots(1, 1, figsize = (15, 10))    # Dipole tilt
    fig_hp   , ax_hp    = plt.subplots(1, 1, figsize = (15, 10))    # Helix-dipole angle
    fig_pmom , ax_pmom  = plt.subplots(1, 1, figsize = (15, 10))    # Electric dipole moment
    
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
            times = times/1000.0
            flat_times = times.to_numpy().flatten()
    
            labelname = simulation_legends[isim] + f" N{iseed}"
            seedoutname = outname + f"_N{iseed}"
    
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
            pdip_moments            = np.zeros(ntimes)
            pdipole_arr             = master_time_df[['p_dipole_x', 'p_dipole_y', 'p_dipole_z']].to_numpy()
            helix_arr               = master_time_df[['helix_global_axis_x', 'helix_global_axis_y', 'helix_global_axis_z']].to_numpy()
            z_axis                  = np.array([0.0, 0.0, 1.0])
            for itime in range(ntimes):
                pdipole_angle = np.arccos(np.dot(pdipole_arr[itime,:], z_axis) / (np.linalg.norm(pdipole_arr[itime,:]) * np.linalg.norm(z_axis)))
                helixpd_angle = np.arccos(np.dot(pdipole_arr[itime,:], helix_arr[itime,:]) / (np.linalg.norm(pdipole_arr[itime,:]) * np.linalg.norm(helix_arr[itime,:])))
                pdipole_tilt_angles[itime] = 180.0 / np.pi * pdipole_angle
                helix_pdipole_angles[itime] = 180.0 / np.pi * helixpd_angle
                pdip_moments[itime] = np.linalg.norm(pdipole_arr[itime,:])
            ax_pdip.plot(times[::stride], pdipole_tilt_angles[::stride], label = labelname)
            ax_hp.plot(times[::stride], helix_pdipole_angles[::stride], label = labelname)
            ax_pmom.plot(times[::stride], pdip_moments[::stride], label = labelname)

            # Get the final values of some things and record them to process later
            # Do some fancy boolean indexing on the times
            final_indices = np.where(flat_times < end_times)
            final_data_dict['zpos'][simtorun]   = np.concatenate((final_data_dict['zpos'][simtorun],    z_pos[final_indices]), axis = 0)
            final_data_dict['helix'][simtorun]  = np.concatenate((final_data_dict['helix'][simtorun],   helicity[final_indices]), axis = 0)
            final_data_dict['tilt'][simtorun]   = np.concatenate((final_data_dict['tilt'][simtorun],    global_tilt[final_indices]), axis = 0)
            final_data_dict['pdip'][simtorun]   = np.concatenate((final_data_dict['pdip'][simtorun],    pdipole_tilt_angles[final_indices]), axis = 0)
            final_data_dict['pmom'][simtorun]   = np.concatenate((final_data_dict['pmom'][simtorun],    pdip_moments[final_indices]), axis = 0)
            final_data_dict['hp'][simtorun]     = np.concatenate((final_data_dict['hp'][simtorun],      helix_pdipole_angles[final_indices]), axis = 0)

            # Plot the initial and final 3d version of the helix director and dipole moment
            # Initial orientations
            fig3dinitial = plt.figure()
            ax3dinitial = fig3dinitial.add_subplot(projection = '3d')
            initial_pdipole_norm = np.linalg.norm(pdipole_arr[0,:])
            ax3dinitial.quiver(0.0, 0.0, 0.0,
                        helix_arr[0,0], helix_arr[0,1], helix_arr[0,2],
                        pivot = 'middle', color = 'tab:blue')
            ax3dinitial.quiver(0.0, 0.0, 0.0,
                        pdipole_arr[0,0]/initial_pdipole_norm, pdipole_arr[0,1]/initial_pdipole_norm, pdipole_arr[0,2]/initial_pdipole_norm,
                        pivot = 'middle', color = 'tab:orange')
            ax3dinitial.set_xlim([-1.0, 1.0])
            ax3dinitial.set_ylim([-1.0, 1.0])
            ax3dinitial.set_zlim([-1.0, 1.0])
            ax3dinitial.set_xlabel("X")
            ax3dinitial.set_ylabel("Y")
            ax3dinitial.set_zlabel("Z")
            fig3dinitial.tight_layout()
            fig3dinitial.savefig('initial/gromacs_initialquiver_' + seedoutname + '_3d.pdf', dpi = fig3dinitial.dpi)
            plt.close()

            fig2dinitial, ax2dinitial = plt.subplots(1, 2)
            ax2dinitial[0].quiver(0.0, 0.0,
                                  helix_arr[0,0], helix_arr[0,1],
                                  units = 'xy', angles = 'xy', scale_units = 'xy', scale = 1, width = 0.05,
                                  pivot = 'middle', color = 'tab:blue')
            ax2dinitial[1].quiver(0.0, 0.0,
                                  helix_arr[0,1], helix_arr[0,2],
                                  units = 'xy', angles = 'xy', scale_units = 'xy', scale = 1, width = 0.05,
                                  pivot = 'middle', color = 'tab:blue')
            ax2dinitial[0].quiver(0.0, 0.0,
                                  pdipole_arr[0,0]/initial_pdipole_norm, pdipole_arr[0,1]/initial_pdipole_norm,
                                  units = 'xy', angles = 'xy', scale_units = 'xy', scale = 1, width = 0.05,
                                  pivot = 'middle', color = 'tab:orange')
            ax2dinitial[1].quiver(0.0, 0.0,
                                  pdipole_arr[0,1]/initial_pdipole_norm, pdipole_arr[0,2]/initial_pdipole_norm,
                                  units = 'xy', angles = 'xy', scale_units = 'xy', scale = 1, width = 0.05,
                                  pivot = 'middle', color = 'tab:orange')
            ax2dinitial[0].set_xlim([-1.0, 1.0])
            ax2dinitial[0].set_ylim([-1.0, 1.0])
            ax2dinitial[0].set_xlabel("X")
            ax2dinitial[0].set_ylabel("Y")
            ax2dinitial[0].set_aspect('equal')
            ax2dinitial[1].set_xlim([-1.0, 1.0])
            ax2dinitial[1].set_ylim([-1.0, 1.0])
            ax2dinitial[1].set_xlabel("Y")
            ax2dinitial[1].set_ylabel("Z")
            ax2dinitial[1].set_aspect('equal')
            fig2dinitial.tight_layout()
            fig2dinitial.savefig('initial/gromacs_initialquiver_' + seedoutname + '_2d.pdf', dpi = fig2dinitial.dpi)
            plt.close()
                    
            # Final orientations
            fig3dfinal = plt.figure()
            ax3dfinal = fig3dfinal.add_subplot(projection = '3d')
            final_pdipole_norm = np.linalg.norm(pdipole_arr[-1,:])
            ax3dfinal.quiver(0.0, 0.0, 0.0,
                        helix_arr[-1,0], helix_arr[-1,1], helix_arr[-1,2],
                        pivot = 'middle', color = 'tab:blue')
            ax3dfinal.quiver(0.0, 0.0, 0.0,
                        pdipole_arr[-1,0]/final_pdipole_norm, pdipole_arr[-1,1]/final_pdipole_norm, pdipole_arr[-1,2]/final_pdipole_norm,
                        pivot = 'middle', color = 'tab:orange')
            ax3dfinal.set_xlim([-1.0, 1.0])
            ax3dfinal.set_ylim([-1.0, 1.0])
            ax3dfinal.set_zlim([-1.0, 1.0])
            ax3dfinal.set_xlabel("X")
            ax3dfinal.set_ylabel("Y")
            ax3dfinal.set_zlabel("Z")
            fig3dfinal.tight_layout()
            fig3dfinal.savefig('final/gromacs_finalquiver_' + seedoutname + '_3d.pdf', dpi = fig3dfinal.dpi)
            plt.close()

            fig2dfinal, ax2dfinal = plt.subplots(1, 2)
            ax2dfinal[0].quiver(0.0, 0.0,
                                helix_arr[-1,0], helix_arr[-1,1],
                                units = 'xy', angles = 'xy', scale_units = 'xy', scale = 1, width = 0.05,
                                pivot = 'middle', color = 'tab:blue')
            ax2dfinal[1].quiver(0.0, 0.0,
                                helix_arr[-1,1], helix_arr[-1,2],
                                units = 'xy', angles = 'xy', scale_units = 'xy', scale = 1, width = 0.05,
                                pivot = 'middle', color = 'tab:blue')
            ax2dfinal[0].quiver(0.0, 0.0,
                                pdipole_arr[-1,0]/final_pdipole_norm, pdipole_arr[-1,1]/final_pdipole_norm,
                                units = 'xy', angles = 'xy', scale_units = 'xy', scale = 1, width = 0.05,
                                pivot = 'middle', color = 'tab:orange')
            ax2dfinal[1].quiver(0.0, 0.0,
                                pdipole_arr[-1,1]/final_pdipole_norm, pdipole_arr[-1,2]/final_pdipole_norm,
                                units = 'xy', angles = 'xy', scale_units = 'xy', scale = 1, width = 0.05,
                                pivot = 'middle', color = 'tab:orange')
            ax2dfinal[0].set_xlim([-1.0, 1.0])
            ax2dfinal[0].set_ylim([-1.0, 1.0])
            ax2dfinal[0].set_xlabel("X")
            ax2dfinal[0].set_ylabel("Y")
            ax2dfinal[0].set_aspect('equal')
            ax2dfinal[1].set_xlim([-1.0, 1.0])
            ax2dfinal[1].set_ylim([-1.0, 1.0])
            ax2dfinal[1].set_xlabel("Y")
            ax2dfinal[1].set_ylabel("Z")
            ax2dfinal[1].set_aspect('equal')
            fig2dfinal.tight_layout()
            fig2dfinal.savefig('final/gromacs_finalquiver_' + seedoutname + '_2d.pdf', dpi = fig2dfinal.dpi)
            plt.close()
                    
            # Categorize by the final global tilt, see what kind of categories we can come up with
            final_tilt_seed = np.mean(global_tilt[final_indices])
            #print(f"  final tilt: {final_tilt_seed}")
            if final_tilt_seed < tilt_criteria or final_tilt_seed > (180.0 - tilt_criteria):
                ax_zpos_tilt.plot(times[::stride], z_pos[::stride], color = 'm')
                ax_helix_tilt.plot(times[::stride], helicity[::stride], color = 'm')
                ax_tilt_tilt.plot(times[::stride], global_tilt[::stride], color = 'm')
                ax_pdip_tilt.plot(times[::stride], pdipole_tilt_angles[::stride], color = 'm')
                ax_hp_tilt.plot(times[::stride], helix_pdipole_angles[::stride], color = 'm')
                ax_pmom_tilt.plot(times[::stride], pdip_moments[::stride], color = 'm')
            else:
                ax_zpos_tilt.plot(times[::stride], z_pos[::stride], color = 'b', alpha = alpha_plot)
                ax_helix_tilt.plot(times[::stride], helicity[::stride], color = 'b', alpha = alpha_plot)
                ax_tilt_tilt.plot(times[::stride], global_tilt[::stride], color = 'b', alpha = alpha_plot)
                ax_pdip_tilt.plot(times[::stride], pdipole_tilt_angles[::stride], color = 'b', alpha = alpha_plot)
                ax_hp_tilt.plot(times[::stride], helix_pdipole_angles[::stride], color = 'b', alpha = alpha_plot)
                ax_pmom_tilt.plot(times[::stride], pdip_moments[::stride], color = 'b', alpha = alpha_plot)
    
        cidx += 1

    #ax_zdist.set_xlabel('Time (ns)')
    #ax_zdist.set_ylabel(r'Z distance ($\AA$)')
    #fig_zdist.tight_layout()
    #fig_zdist.savefig('gromacs_zdist_' + outname + '.pdf', dpi = fig_zdist.dpi)
    #plt.close()
   
    # zpos
    ax_zpos.set_xlabel('Time (ns)')
    ax_zpos.set_ylabel(r'z ($\AA$)')
    ax_zpos.set_ylim([-30.0, 30.0])
    fig_zpos.tight_layout()
    fig_zpos.savefig('gromacs_zpos_' + outname + '.pdf', dpi = fig_zpos.dpi)
    plt.close()

    # helicity 
    ax_helix.set_xlabel('Time (ns)')
    ax_helix.set_ylabel('Helical nature (AU)')
    ax_helix.set_ylim([0.0, 1.05])
    fig_helix.tight_layout()
    fig_helix.savefig('gromacs_helicity_' + outname + '.pdf', dpi = fig_helix.dpi)
    plt.close()

    # global tilt
    ax_tilt.set_xlabel('Time (ns)')
    ax_tilt.set_ylabel('Global tilt (deg)')
    ax_tilt.set_ylim([0.0, 180.0])
    fig_tilt.tight_layout()
    fig_tilt.savefig('gromacs_tilt_' + outname + '.pdf', dpi = fig_tilt.dpi)
    plt.close()
   
    # helix dipole tilt
    ax_pdip.set_xlabel('Time (ns)')
    ax_pdip.set_ylabel('Helix dipole tilt (deg)')
    ax_pdip.set_ylim([0.0, 180.0])
    fig_pdip.tight_layout()
    fig_pdip.savefig('gromacs_pdipoletilt_' + outname + '.pdf', dpi = fig_pdip.dpi)
    plt.close()
   
    # helix dipole helix angle 
    ax_hp.set_xlabel('Time (ns)')
    ax_hp.set_ylabel('Helix dipole - Helix axis angle (deg)')
    ax_hp.set_ylim([0.0, 180.0])
    fig_hp.tight_layout()
    fig_hp.savefig('gromacs_helixpdipoleangle_' + outname + '.pdf', dpi = fig_hp.dpi)
    plt.close()

    # pdipole moment
    ax_pmom.set_xlabel('Time (ns)')
    ax_pmom.set_ylabel('Helix electric dipole moment (?)')
    ax_pmom.set_ylim([0.0, 50.0])
    fig_pmom.tight_layout()
    fig_pmom.savefig('gromacs_pmoment_' + outname + '.pdf', dpi = fig_pmom.dpi)
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
    print(f"      {yvals} +/- {yerrs}")
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

# Also the categorized tilts and things like that!
ax_zpos_tilt.set_xlabel('Time (ns)')
ax_zpos_tilt.set_ylabel(r'z ($\AA$)')
ax_zpos_tilt.set_ylim([-30.0, 30.0])
fig_zpos_tilt.tight_layout()
fig_zpos_tilt.savefig('gromacs_zpos_tilt.pdf', dpi = fig_zpos_tilt.dpi)
plt.close()

ax_helix_tilt.set_xlabel('Time (ns)')
ax_helix_tilt.set_ylabel('Helical nature (AU)')
ax_helix_tilt.set_ylim([0.0, 1.05])
fig_helix_tilt.tight_layout()
fig_helix_tilt.savefig('gromacs_helicity_tilt.pdf', dpi = fig_helix_tilt.dpi)
plt.close()
    
ax_tilt_tilt.set_xlabel('Time (ns)')
ax_tilt_tilt.set_ylabel('Global tilt (deg)')
ax_tilt_tilt.set_ylim([0.0, 180.0])
fig_tilt_tilt.tight_layout()
fig_tilt_tilt.savefig('gromacs_tilt_tilt.pdf', dpi = fig_tilt_tilt.dpi)
plt.close()
   
ax_pdip_tilt.set_xlabel('Time (ns)')
ax_pdip_tilt.set_ylabel('Helix dipole tilt (deg)')
ax_pdip_tilt.set_ylim([0.0, 180.0])
fig_pdip_tilt.tight_layout()
fig_pdip_tilt.savefig('gromacs_pdipoletilt_tilt.pdf', dpi = fig_pdip_tilt.dpi)
plt.close()
   
ax_hp_tilt.set_xlabel('Time (ns)')
ax_hp_tilt.set_ylabel('Helix dipole - Helix axis angle (deg)')
ax_hp_tilt.set_ylim([0.0, 180.0])
fig_hp_tilt.tight_layout()
fig_hp_tilt.savefig('gromacs_helixpdipoleangle_tilt.pdf', dpi = fig_hp_tilt.dpi)
plt.close()

ax_pmom_tilt.set_xlabel('Time (ns)')
ax_pmom_tilt.set_ylabel('Helix electric dipole moment (?)')
ax_pmom_tilt.set_ylim([0.0, 50.0])
fig_pmom_tilt.tight_layout()
fig_pmom_tilt.savefig('gromacs_pmoment_tilt.pdf', dpi = fig_pmom_tilt.dpi)
plt.close()

