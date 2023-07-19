#!/usr/bin/env python3

# XXX: Put a license here

""" Figure 1 hardcoded oneoff script """

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
from scipy.signal import savgol_filter

# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src', 'lib'))
from stylelib.common_styles import *
from common import create_datadir

# Set the style for the plots
plt.style.use(septin_poster_stl)

# Move thit into a subclass at some point
def shaded_error(ax, x, y, error, alpha, color, linestyle = 'solid', label = None, linewidth = 1):
    r""" Plot a shaded error bar with colors
    """
    ax.plot(x, y, color = color, linestyle = linestyle, label = label, linewidth = linewidth)
    ax.fill_between(x, y-error, y+error,
                    alpha = alpha, edgecolor = color, facecolor = color, linestyle = linestyle, linewidth = linewidth)

# Doesn't work and can't figure out why right now
#def convert_pdf_gs(name, newname):
#    r""" Convert the PDF file using GS to make it smaller and more manageable
#    """
#    pdf_p = subprocess.Popen(["gs", "-sDevice=pdfwrite", "-dCompatibilityLevel=1.4", "-dNOPAUSE", "-dQUIET", "-dBATCH", "-sOutputFile={}".format(newname), name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
#    pdf_output, pdf_errors = pdf_p.communicate()
#    pdf_p.wait()
#    print(pdf_output.decode("utf-8"))
#    print(pdf_errors.decode("utf-8"))

#def gs_opt(filename):
#    filenameTmp = filename.split('.')[-2]+'_tmp.pdf'
#    gs = ['gs',
#          '-sDEVICE=pdfwrite',
#          '-dEmbedAllFonts=false',
#          '-dSubsetFonts=true',             # Create font subsets (default)
#          '-dPDFSETTINGS=/prepress',        # Image resolution
#          '-dDetectDuplicateImages=true',   # Embeds images used multiple times only once
#          '-dCompressFonts=true',           # Compress fonts in the output (default)
#          '-dNOPAUSE',                      # No pause after each image
#          '-dQUIET',                        # Suppress output
#          '-dBATCH',                        # Automatically exit
#          '-sOutputFile='+filenameTmp,      # Save to temporary output
#          filename]                         # Input file
#    
#    subprocess.run(gs)                                      # Create temporary file
#    #subprocess.run(['rm', filename], shell=True)            # Delete input file
#    #subprocess.run(['mv', filenameTmp, filename], shell=True) # Rename temporary to input file
#    os.rename(filenameTmp, filename)


main_path = os.path.abspath('/Users/cedelmaier/Projects/Biophysics/septin_project/atomistic/simulations/data/')
external_path = os.path.abspath('/Volumes/T7/data/septin_project/datasets/')

# Create some external data directories
allseed_datadir = create_datadir(os.getcwd(), datadir_name = "allseeds")
mean_datadir = create_datadir(os.getcwd(), datadir_name = "mean")

#simnames = {
#            "neutral_fold_00": "rfmonomer_aglipid_11x11_zdepth00_rotx0_50mMKCl_long",
#            "charge_fold_00": "unbiased_zdepth00_rotx0_helix_50mMKCl",
#            "charge_fold_15": "agmonomer_aglipid_11x11_zdepth15_rotx0_50mMKCl",
#            "charge_fold_30": "unbiased_50mMKCl_solution",
#            "charge_unfold_15": "agmonomermelt_aglipid_11x11_zdepth15_rotx0_50mMKCl",
#            }
simnames = {
            "neutral_fold_00":      os.path.join(main_path, "rfmonomer_aglipid_11x11_zdepth00_rotx0_50mMKCl_long"),
            "neutral_fold_15":      os.path.join(external_path, "rfmonomer_aglipid_11x11_zdepth15_50mMKCl"),
            "neutral_fold_30":      os.path.join(external_path, "rfmonomer_aglipid_11x11_zdepth30_50mMKCl"),
            "neutral_unfold_15":    os.path.join(external_path, "rfmonomermelt_aglipid_11x11_zdepth15_rotx0_50mMCKl"),
            "charge_fold_00":       os.path.join(external_path, "gromacs_zdepth00_rotx0_helix_50mMKCl"),
            "charge_fold_15":       os.path.join(main_path, "agmonomer_aglipid_11x11_zdepth15_rotx0_50mMKCl"),
            "charge_fold_30":       os.path.join(external_path, "unbiased_50mMKCl_solution"),
            "charge_unfold_15":     os.path.join(main_path, "agmonomermelt_aglipid_11x11_zdepth15_rotx0_50mMKCl"),
            }

#seednames = ["N1"]
seednames = ["N1", "N2", "N3", "N4"]

axis_names = {}
#axis_names['zpos']  = r"z ($\AA$)"
#axis_names['helix'] = r"Fraction helix (AU)"
axis_names['zpos'] = r"z (nm)"
axis_names['alpharmsd'] = r"Helical content (AU)"

# Load up the ylow ad yhi dicts
ylow_dict = {"zpos": 0.0,
             "alpharmsd": 0.0,
             "helix": 0.0,
             "tilt": 0.0,
             "pdip": 0.0,
             "pmom": 0.0,
             "hp": 0.0,
             "zforce": -8000.0,
             "perptorque": -800000.0}
yhi_dict = {"zpos": 10.0,
            "alpharmsd": 13.1,
            "helix": 1.05,
            "tilt": 180.0,
            "pdip": 180.0,
            "pmom": 40.0,
            "hp": 180.0,
            "zforce": 8000.0,
            "perptorque": 800000.0}
# Colorblind color cycle
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

# Custom plot options
plt.style.use(septin_paper_stl)

# Colorblind index
cidx = 0
for simname,seedname in simnames.items():
    # print the simulation we are doing
    print(simname, seedname)

    fig_zpos, axarr_zpos            = plt.subplots(1, 2, figsize = (5.0, 2.5), sharey = True)
    fig_zpos_mean, axarr_zpos_mean  = plt.subplots(1, 2, figsize = (5.0, 2.5), sharey = True)
    fig_heli, axarr_heli            = plt.subplots(1, 2, figsize = (5.0, 2.5), sharey = True)
    fig_heli_mean, axarr_heli_mean  = plt.subplots(1, 2, figsize = (5.0, 2.5), sharey = True)

    # pretend we are a seed
    cidx = 0
    leaf0_list = []
    leaf1_list = []
    leaf1_list_top = []
    center_list = []
    heli_list = []
    for sd in seednames:
        file_path_trajectory = os.path.join(seedname, sd, sd + ".h5")
        master_time_df = pd.read_hdf(file_path_trajectory)
        print(master_time_df)

        # Z position
        # Conver the Z positions to nm rather than angstroms
        z_protein   = master_time_df[['helix_z']].to_numpy().flatten()/10.0
        z_leaf0     = master_time_df[['leaflet0_z']].to_numpy().flatten()/10.0
        z_leaf1     = master_time_df[['leaflet1_z']].to_numpy().flatten()/10.0
        z_lipid     = master_time_df[['lipid_z']].to_numpy().flatten()/10.0
        z_pbc       = master_time_df[['unit_cell_z']].to_numpy().flatten()/10.0
        # Subtract off the position of the lipid COM from everybody else
        pbc_full_z = z_pbc
        z_protein = z_protein - z_lipid
        z_leaf0 = z_leaf0 - z_lipid
        z_leaf1 = z_leaf1 - z_lipid
        z_pbc = z_pbc/2 - z_lipid
        # Correc the position if under the lower leaflet
        for idx in range(len(z_protein)):
            if z_protein[idx] < z_leaf1[idx]:
                z_protein[idx] = z_pbc[idx] - z_protein[idx]

        xdata = master_time_df.index / 1000.0;
        ydata = z_protein

        max_time = xdata[-1]

        # This version puts the breaks in the x-axis
        plt.figure(fig_zpos)
        axarr_zpos[0].plot(xdata, ydata, linewidth = 1, color = CB_color_cycle[cidx])
        axarr_zpos[1].plot(xdata, ydata, linewidth = 1, color = CB_color_cycle[cidx])

        leaf0_list.append(z_leaf0)
        leaf1_list.append(z_leaf1)
        leaf1_list_top.append(z_leaf1 + pbc_full_z)
        center_list.append(ydata)

        # Do the helix as well
        plt.figure(fig_heli)
        #helicity    = master_time_df[['helicity']].to_numpy().flatten() 
        helicity    = master_time_df[['plumed_alpha']].to_numpy().flatten() 
        axarr_heli[0].plot(xdata, helicity, linewidth = 1, color = CB_color_cycle[cidx])
        axarr_heli[1].plot(xdata, helicity, linewidth = 1, color = CB_color_cycle[cidx])
        heli_list.append(helicity)

        cidx += 1

    # Do the leaflet limits
    plt.figure(fig_zpos)
    leaf0_mean = np.mean(np.array(leaf0_list), axis=0)
    leaf1_mean = np.mean(np.array(leaf1_list), axis=0)
    leaf1_top_mean = np.mean(np.array(leaf1_list_top), axis=0)
    leaf0_std = np.std(np.array(leaf0_list), axis=0, ddof=1)
    leaf1_std = np.std(np.array(leaf1_list), axis=0, ddof=1)
    leaf1_top_std = np.std(np.array(leaf1_list_top), axis=0, ddof=1)

    # Mean version as well
    center_mean = np.mean(np.array(center_list), axis=0)
    center_std  = np.std(np.array(center_list), axis=0, ddof=1)
    heli_mean   = np.mean(np.array(heli_list), axis=0)
    heli_std    = np.std(np.array(heli_list), axis=0, ddof=1)

    shaded_error(axarr_zpos[0], xdata, leaf0_mean, leaf0_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos[1], xdata, leaf0_mean, leaf0_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos[0], xdata, leaf1_mean, leaf1_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos[1], xdata, leaf1_mean, leaf1_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos[0], xdata, leaf1_top_mean, leaf1_top_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos[1], xdata, leaf1_top_mean, leaf1_top_std, alpha = 0.5, color = 'slategrey')

    # Mean version duplicates lots of work
    shaded_error(axarr_zpos_mean[0], xdata, leaf0_mean, leaf0_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos_mean[1], xdata, leaf0_mean, leaf0_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos_mean[0], xdata, leaf1_mean, leaf1_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos_mean[1], xdata, leaf1_mean, leaf1_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos_mean[0], xdata, leaf1_top_mean, leaf1_top_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos_mean[1], xdata, leaf1_top_mean, leaf1_top_std, alpha = 0.5, color = 'slategrey')
    shaded_error(axarr_zpos_mean[0], xdata, center_mean, center_std, alpha = 0.5, color = CB_color_cycle[0])
    shaded_error(axarr_zpos_mean[1], xdata, center_mean, center_std, alpha = 0.5, color = CB_color_cycle[0])

    # Set the ylimits
    axarr_zpos[0].set_ylim(ylow_dict["zpos"], yhi_dict["zpos"])
    axarr_zpos[1].set_ylim(ylow_dict["zpos"], yhi_dict["zpos"])
    axarr_zpos_mean[0].set_ylim(ylow_dict["zpos"], yhi_dict["zpos"])
    axarr_zpos_mean[1].set_ylim(ylow_dict["zpos"], yhi_dict["zpos"])

    # Set the limits of the two subplots differently
    axarr_zpos[0].set_xlim(0.0, 10.0)
    axarr_zpos[1].set_xlim(max_time - 100.0, max_time)
    axarr_zpos_mean[0].set_xlim(0.0, 10.0)
    axarr_zpos_mean[1].set_xlim(max_time - 100.0, max_time)

    # Set so we can't see the spines, etc
    axarr_zpos[0].spines['right'].set_visible(False)
    axarr_zpos[0].tick_params(labelright = False)
    axarr_zpos_mean[0].spines['right'].set_visible(False)
    axarr_zpos_mean[0].tick_params(labelright = False)

    axarr_zpos[1].spines['left'].set_visible(False)
    axarr_zpos[1].tick_params(labelleft = False)
    axarr_zpos[1].yaxis.tick_right()
    axarr_zpos[1].tick_params(right=False)
    axarr_zpos_mean[1].spines['left'].set_visible(False)
    axarr_zpos_mean[1].tick_params(labelleft = False)
    axarr_zpos_mean[1].yaxis.tick_right()
    axarr_zpos_mean[1].tick_params(right=False)

    # Create the breaks
    d = 0.01
    kwargs = dict(transform=axarr_zpos[0].transAxes, color = 'k', clip_on = False)
    axarr_zpos[0].plot((1-d, 1+d), (-d, d), linewidth=1, **kwargs)
    axarr_zpos[0].plot((1-d, 1+d), (1-d, 1+d), linewidth=1, **kwargs)
    kwargs.update(transform=axarr_zpos[1].transAxes)
    axarr_zpos[1].plot((-d, d), (-d, d), linewidth=1, **kwargs)
    axarr_zpos[1].plot((-d, d), (1-d, 1+d), linewidth=1, **kwargs)
    kwargs = dict(transform=axarr_zpos_mean[0].transAxes, color = 'k', clip_on = False)
    axarr_zpos_mean[0].plot((1-d, 1+d), (-d, d), linewidth=1, **kwargs)
    axarr_zpos_mean[0].plot((1-d, 1+d), (1-d, 1+d), linewidth=1, **kwargs)
    kwargs.update(transform=axarr_zpos_mean[1].transAxes)
    axarr_zpos_mean[1].plot((-d, d), (-d, d), linewidth=1, **kwargs)
    axarr_zpos_mean[1].plot((-d, d), (1-d, 1+d), linewidth=1, **kwargs)

    spaceshift = 0.10
    fig_zpos.supxlabel(r"Time (ns)", x = 0.5+spaceshift/2.0)
    axarr_zpos[0].set_ylabel(axis_names["zpos"])
    fig_zpos_mean.supxlabel(r"Time (ns)", x = 0.5+spaceshift/2.0)
    axarr_zpos_mean[0].set_ylabel(axis_names["zpos"])

    # Layout options must come in order
    plt.figure(fig_zpos)
    fig_zpos.tight_layout()
    plt.subplots_adjust(wspace=spaceshift, bottom=0.16)
    allseed_zpos_name = allseed_datadir + "/figure1_zpos_" + simname + "_allseeds.png"
    plt.savefig(allseed_zpos_name, dpi = 600)
    plt.figure(fig_zpos_mean)
    fig_zpos_mean.tight_layout()
    plt.subplots_adjust(wspace=spaceshift, bottom=0.16)
    mean_zpos_name = mean_datadir + "/figure1_zpos_" + simname + "_mean.png"
    plt.savefig(mean_zpos_name, dpi = 600)

    # Do the helicity too
    shaded_error(axarr_heli_mean[0], xdata, heli_mean, heli_std, alpha = 0.5, color = CB_color_cycle[0])
    shaded_error(axarr_heli_mean[1], xdata, heli_mean, heli_std, alpha = 0.5, color = CB_color_cycle[0])
    plt.figure(fig_heli)
    # Set the ylimits
    axarr_heli[0].set_ylim(ylow_dict["alpharmsd"], yhi_dict["alpharmsd"])
    axarr_heli[1].set_ylim(ylow_dict["alpharmsd"], yhi_dict["alpharmsd"])
    axarr_heli_mean[0].set_ylim(ylow_dict["alpharmsd"], yhi_dict["alpharmsd"])
    axarr_heli_mean[1].set_ylim(ylow_dict["alpharmsd"], yhi_dict["alpharmsd"])

    # Set the limits of the two subplots differently
    axarr_heli[0].set_xlim(0.0, 10.0)
    axarr_heli[1].set_xlim(max_time - 100.0, max_time)
    axarr_heli_mean[0].set_xlim(0.0, 10.0)
    axarr_heli_mean[1].set_xlim(max_time - 100.0, max_time)

    # Set so we can't see the spines, etc
    axarr_heli[0].spines['right'].set_visible(False)
    axarr_heli[0].tick_params(labelright = False)
    axarr_heli_mean[0].spines['right'].set_visible(False)
    axarr_heli_mean[0].tick_params(labelright = False)

    axarr_heli[1].spines['left'].set_visible(False)
    axarr_heli[1].tick_params(labelleft = False)
    axarr_heli[1].yaxis.tick_right()
    axarr_heli[1].tick_params(right=False)
    axarr_heli_mean[1].spines['left'].set_visible(False)
    axarr_heli_mean[1].tick_params(labelleft = False)
    axarr_heli_mean[1].yaxis.tick_right()
    axarr_heli_mean[1].tick_params(right=False)

    # Create the breaks
    d = 0.01
    kwargs = dict(transform=axarr_heli[0].transAxes, color = 'k', clip_on = False)
    axarr_heli[0].plot((1-d, 1+d), (-d, d), linewidth=1, **kwargs)
    axarr_heli[0].plot((1-d, 1+d), (1-d, 1+d), linewidth=1, **kwargs)
    kwargs.update(transform=axarr_heli[1].transAxes)
    axarr_heli[1].plot((-d, d), (-d, d), linewidth=1, **kwargs)
    axarr_heli[1].plot((-d, d), (1-d, 1+d), linewidth=1, **kwargs)
    kwargs = dict(transform=axarr_heli_mean[0].transAxes, color = 'k', clip_on = False)
    axarr_heli_mean[0].plot((1-d, 1+d), (-d, d), linewidth=1, **kwargs)
    axarr_heli_mean[0].plot((1-d, 1+d), (1-d, 1+d), linewidth=1, **kwargs)
    kwargs.update(transform=axarr_heli_mean[1].transAxes)
    axarr_heli_mean[1].plot((-d, d), (-d, d), linewidth=1, **kwargs)
    axarr_heli_mean[1].plot((-d, d), (1-d, 1+d), linewidth=1, **kwargs)

    fig_heli.supxlabel(r"Time (ns)", x = 0.5+spaceshift/2.0)
    axarr_heli[0].set_ylabel(axis_names["alpharmsd"])
    fig_heli_mean.supxlabel(r"Time (ns)", x = 0.5+spaceshift/2.0)
    axarr_heli_mean[0].set_ylabel(axis_names["alpharmsd"])

    # Layout options must come in order
    plt.figure(fig_heli)
    fig_heli.tight_layout()
    plt.subplots_adjust(wspace=spaceshift, bottom=0.16)
    allseed_heli_name = allseed_datadir + "/figure1_alpharmsd_" + simname + "_allseeds.png"
    plt.savefig(allseed_heli_name, dpi = 600)

    plt.figure(fig_heli_mean)
    fig_heli_mean.tight_layout()
    plt.subplots_adjust(wspace=spaceshift, bottom=0.16)
    mean_heli_name = mean_datadir + "/figure1_alpharmsd_" + simname + "_mean.png"
    plt.savefig(mean_heli_name, dpi = 600)

    # This is a trick to get the PDF files of the right size, at MPL decides to make them rather
    # unpleasant to use, in terms of size, for things like Illustrator, Keynote, etc.
    #convert_pdf_gs(allseed_zpos_name, allseed_zpos_name_gs)
    #convert_pdf_gs(mean_zpos_name, mean_zpos_name_gs)
    #convert_pdf_gs(allseed_heli_name, allseed_heli_name_gs)
    #convert_pdf_gs(mean_heli_name, mean_heli_name_gs)

    plt.close('all')

    # For each simname, save off the zpos and alpha that we've come up with in 
    # pandas format so that we don't have to go through this repeatedly for 
    # actually plotting the thing
    simdict = {
            "timepoints": xdata,
            "leaf0_mean": leaf0_mean,
            "leaf0_std": leaf0_std,
            "leaf1_mean": leaf1_mean,
            "leaf1_std": leaf1_std,
            "leaf1_top_mean": leaf1_top_mean,
            "leaf1_top_std": leaf1_top_std,
            "center_mean": center_mean,
            "center_std": center_std,
            "heli_mean": heli_mean,
            "heli_std": heli_std
            }
    print(simdict)
    simdf = pd.DataFrame(simdict)
    print(simdf)
    simdf.to_hdf(os.path.join("./", "{}.h5".format(simname)), key = "simname", mode = "w")
