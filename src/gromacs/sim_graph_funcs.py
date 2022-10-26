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
from common import *

# General histogram functions
def graph_sim_histogram(xdata,
                        ax,
                        mtitle = '',
                        xtitle = '',
                        ytitle = '',
                        xrange = (0.0, 100.0),
                        color = 'b'):
    r""" General histogram plotting tool, with errorbars
    """
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

    # Get the actual original histogram
    hist, bin_edges = np.histogram(xdata, range = xrange)
    bin_mids = moving_average(bin_edges)
    n_points = xdata.size

    #ax.hist(xdata, range=xrange, color = color)
    ax.scatter(x = bin_mids, y = hist, color = color, marker = "s")

    return (hist, bin_mids, n_points)

# Scatter plot bewtween results
def graph_sim_scatter(xdata,
                      ydata,
                      ax,
                      mtitle = '',
                      xtitle = '',
                      ytitle = '',
                      color = 'b'):
    r""" General scatter plot for two variables
    """
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

    ax.scatter(xdata, ydata, color = color)
