# Something

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def graph_scatter(sdlabel,
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

