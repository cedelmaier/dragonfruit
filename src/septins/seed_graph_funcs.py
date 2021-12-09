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
    #ax.set_title(mtitle)
    #ax.set_ylabel(ytitle)
    #if xlabel: ax.set_xlabel(xtitle)

    # XXX Max frame hardcoded to 200
    #hist, bin_edges = np.histogram(df['free'], bins = 10, range = (0, 200))
    #hist = df.hist(ax = axarr, bins = 20, range = (0, 200))
    #axarr.set_xlabel("Lifetime (frame)")

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
