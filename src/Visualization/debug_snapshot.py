#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 16:46:21 2017

@author: duman
"""

""" Plot snapshots of timeframes from the simulation"""
# NOTE: This script only works for filaments at the moment.

##############################################################################

import sys
sys.path.append('../../include')

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import h5py
import os
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import matplotlib.colors as mplcolors

import read_write
import misc_tools
import data_separator
import plotter
import data_structures

import seaborn as sns
sns.set(style="white",context='paper',
        font_scale=1.2,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 30,
            'font.family': 'sans',"figure.dpi":300,
            "xtick.major.size": 8, "ytick.major.size": 8,
            'grid.linestyle': '--'})

##############################################################################

def get_args():
    """ get the command line arguments"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-fd", "--folder", \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Bidisperse_Filaments/.../")
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Bidisperse_Filaments/IMAGES/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Bidisperse_Filaments/IMAGES/")
    parser.add_argument("-ti","--init_time", nargs="?", const=100, type=int, \
                        help="First frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=1900, type=int, \
                        help="Last frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-c","--colorid", type=str, \
                        help ="Decide on the coloring -id or orient-")
    parser.add_argument("-st", "--simtype", \
                        help="Simulation type --filaments or cells--")
    parser.add_argument("-s","--savepdf", action="store_true", \
                        help ="Decide whether to save as pdf or not")
    args = parser.parse_args()

    return args

##############################################################################

def set_color_pallette(beads, sim, colorid):
    """ set the color pallette"""

    quant_steps = 2056
    if colorid == "id":
        minval = 0.
        maxval = sim.npols
        cmap_ax = plt.cm.get_cmap('jet', quant_steps)
        norm_ax = mpl.colors.Normalize(vmin=minval, vmax=maxval)
        cidx = beads.pid
    elif colorid == "orient":
        minval = 0.
        maxval = 2.*np.pi
        cmap_ax = plt.cm.get_cmap('hsv', quant_steps)
        norm_ax = mpl.colors.Normalize(vmin=minval, vmax=maxval)
        cidx = beads.ori

    return cidx, minval, maxval, cmap_ax, norm_ax

##############################################################################

def set_plot_props(sim):
    """ set plot properties"""

    full_box_downlim = -1.
    full_box_uplim = sim.lx+1.

    return full_box_downlim/sim.bond_length, full_box_uplim/sim.bond_length

##############################################################################

def plot_frames(beads, sim, ti, tf, savebase, colorid, savepdf):
    """ plot the selected timeframes"""

    p = plotter.GeneralPlot()
    ax0 = p.ax0
    savepath = "/usr/users/iff_th2/duman/Desktop/debug.pdf"

    rcom = np.zeros((sim.nsteps), dtype=np.float64)
    rcomi = np.zeros((sim.nsteps), dtype=np.float64)
    time = np.arange(sim.nsteps)*sim.dt

    for step in np.arange(0, sim.nsteps):

        print "step / nsteps : ", str(step), " / ", str(sim.nsteps)

        xcom = np.mean(beads.xu[step,0,:])
        ycom = np.mean(beads.xu[step,1,:])
        rcom[step] = np.sqrt(xcom**2 + ycom**2)
        xcomi = np.mean(beads.xi[step,0,:])
        ycomi = np.mean(beads.xi[step,1,:])
        rcomi[step] = np.sqrt(xcomi**2 + ycomi**2)

    ax0.plot(time, rcom, linewidth=2.0, color="red", label="unwrapped")
    ax0.plot(time, rcomi, '--', linewidth=2.0, color="blue", label="wrapped")
    plt.savefig(savepath, dpi=300, bbox_inches='tight', pad_inches=0.08)

    p.fig.clf()
    plt.close()

    return

##############################################################################

def read_data(folder, simtype):

    beads = data_structures.Beads(folder)
    pols = data_structures.Polymers(folder, simtype)
    if simtype == "filaments":
        sim = data_structures.SimulationFilaments(folder)
    elif simtype == "cells":
        sim = data_structures.SimulationCells(folder)
    else:
        raise ValueError("simtype has to be filaments or cells!")

    return beads, pols, sim

##############################################################################

def main():

    ### get the command line arguments

    args = get_args()

    ### read the data and general information from the folder

    beads, pols, sim = read_data(args.folder, args.simtype)
    print "folder = ", args.folder

    print "Calculating image positions of beads"
    beads.calc_img_positions(sim.lx)

    ### plot the data in the given time window

    plot_frames(beads, sim, args.init_time, args.fin_time, \
                args.savebase, args.colorid, args.savepdf)

    return

##############################################################################

if __name__ == '__main__':
    main()

##############################################################################

