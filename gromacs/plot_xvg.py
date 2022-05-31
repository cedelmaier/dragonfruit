#!/usr/bin/env python3

""" File to plot XVG files since we have issues """

import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np


# Magic to get the library directory working properly
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src', 'lib'))
from stylelib.common_styles import septin_runs_stl

xvgname = sys.argv[1]

mtitle = ""
xlabel = ""
ylabel = ""

x, y = [], []

with open(xvgname) as f:
    for line in f:
        if line.startswith("#"):
            continue
        elif line.startswith("@"):
            cols = line.split()

            if cols[1] == "title":
                mtitle = line.split("title")[1].strip()
            elif cols[1] == "xaxis":
                xlabel = line.split("label")[1].strip()
            elif cols[1] == "yaxis":
                ylabel = line.split("label")[1].strip()
        elif re.search("^[0-9]", line):
            # We have what we want
            x.append(np.float32(line.split()[0]))
            y.append(np.float32(line.split()[1]))

x = np.array(x)
y = np.array(y)

# Clean up quotation marks
mtitle = mtitle[1:-1]
xlabel = xlabel[1:-1]
ylabel = ylabel[1:-1]

plt.style.use(septin_runs_stl)
fig, ax = plt.subplots(1,1, figsize=(15,10))
ax.set_title(mtitle)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.plot(x, y, color = "k")
fig.tight_layout()

savename = xvgname.split('.')[0] + '.pdf'
fig.savefig(savename, dpi=fig.dpi)
