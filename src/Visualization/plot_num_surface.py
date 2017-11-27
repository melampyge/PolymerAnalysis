#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 17:05:32 2017

@author: duman
"""

""" Plot number of surface cells as a function of a simulation parameter"""

###

##############################################################################

import sys
sys.path.append('../../include')

import numpy as np
import argparse
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import plotter
import read_write
import misc_tools
import data_separator

##############################################################################

def get_args(read_fnc):
    """ get the command line arguments"""

    parser = argparse.ArgumentParser()

    ### type info

    parser.add_argument("-st", "--sim_type", type=str, \
                        help="Simulation type --cells, bifils, highfils--")
    parser.add_argument("-pt", "--plot_type", type=str, \
                        help="Plot type --2D, 1D, phase, single--")
    parser.add_argument("-an", "--analysis_name", type=str, \
                        help="Analysis name (also, folder+file names) --MSD_subt--")

    ### param info

    parser.add_argument("-lg", "--legend", type=str, nargs="?", \
                        help="Legend parameter")
    parser.add_argument("-ct", "--control", type=str, nargs="?", \
                        help="Control parameter(s)")

    ### param data

    parser.add_argument("-d", "--density", type=float, nargs="?", \
                        help="System density in packing fraction")
    parser.add_argument("-k", "--kappa", type=float, nargs="?", \
                        help="Bending rigidity")
    parser.add_argument("-f", "--fp", type=float, nargs="?", \
                        help="Propulsion force")
    parser.add_argument("-e", "--eps", type=float, nargs="?", \
                        help="Strength of LJ potential")
    parser.add_argument("-a", "--areak", type=float, nargs="?", \
                        help="Area constraint potential strength")

    ### general info

    parser.add_argument("-s","--savepdf", action="store_true", \
                        help="Decide whether to save in pdf or not")
    args = parser.parse_args()

    separator = data_separator.Separator(args, read_fnc)

    return separator, args

##############################################################################

def power_law(x, a, b):

    return a * x**b

##############################################################################

class PlotProps:
    """ encapsulate plot properties"""

    def __init__(self):

        self.read_fnc = read_write.read_single_analysis_data
        self.num_ticks = 5
        self.ax_len = 1.0                          # Length of one subplot square box
        self.ax_b = 0.0                            # Beginning/offset of the subplot in the box
        self.ax_sep = 0.0                          # Separation length between two subplots
        self.total_subplots_in_x = 1               # Total number of subplots

        return

##############################################################################

def restructure_keys(keys):
    """ restructure double keys into single keys"""

    lg_keys = []
    ct_keys = []
    for key in keys:
        if key[0] not in lg_keys:
            lg_keys.append(key[0])
        if key[1] not in ct_keys:
            ct_keys.append(key[1])

    lg_keys = np.sort(lg_keys)
    ct_keys = np.sort(ct_keys)

    return lg_keys, ct_keys

##############################################################################

def select_data(data, lg_keys, ct_keys):
    """ select and/or normalize the data for plotting"""

    xd = {}
    yd = {}
    simd = {}
    for j, lkey in enumerate(lg_keys):
        xl = []
        yl = []
        sims = []
        for i, ckey in enumerate(ct_keys):
            key = (lkey, ckey)
            sims.append(data[key].sim)
            if (data[key].sim.pe != -1):
                yl.append(data[key].data)
                xl.append(data[key].sim.pe)
        xd[lkey] = np.array(xl)
        for j in range(len(yl)):
            if yl[j] < 1e-30:
                yl[j] = 0.0
        yd[lkey] = np.array(yl)
        simd[lkey] = sims

    return xd, yd, simd, lg_keys

##############################################################################

@plotter.plot_1d
def plot_analysis(data, sep, savepdf, *args):
    """ plot the analysis data specific for this analysis"""

    fig = args[0]       # figure handle
    ax0 = args[1]       # axes handle

    ### normalize or select data

    keys = data.keys()
    lg_keys, ct_keys = restructure_keys(keys)
    xd, yd, simd, keys = select_data(data, lg_keys, ct_keys)
    gsim = simd[simd.keys()[0]][0]
    avg_radius = gsim.avg_radius
    data_separator.assign_physicalvalues(sep.fixed_param, gsim)

    ### plot the data


    for j, key in enumerate(keys):

        x = np.array(xd[key])
        y = np.array(yd[key])
        sim = simd[key]

        sep.legend_param.assign_physicalvalue(sim[0])
        label = data_separator.gen_label(sep.legend_param)
        print sim[0].eps, sim[0].pe, y
        line0 = ax0.plot(x, y, \
                        label=label, linewidth=3.0)
        #line1 = ax0.scatter(x, y)

    ### scales

    ax0.set_xscale('log')
    ax0.set_yscale('log')

    ### title

    title = data_separator.gen_title(sep.fixed_param)
    ax0.set_title(title, fontsize=20)

    ### side labels

    ax0.set_xlabel(r'Pe', fontsize=30)
    ax0.set_ylabel(r'$\mathrm{N}_{\mathrm{surface}}$', fontsize=30)

    ### limits

    #ax0.set_xlim((0, gsim.lx/2./avg_radius))
    #ax0.set_ylim(())

    ### ticks

    #ax0.xaxis.set_ticks(np.linspace(0, 15, num_ticks, endpoint=True))
    #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    ax0.tick_params(axis='both', which='major', labelsize=20)

    ### legend

    ax0.legend(bbox_to_anchor=(0.05, 0., 0.65, 1.), loc=1, borderaxespad=0., \
        prop={'size': 20}, mode="expand", frameon=False)

    return fig

##############################################################################

def main():

    props = PlotProps()
    separator, args = get_args(props.read_fnc)
    plot_analysis(separator.data, separator, args.savepdf)

    return

##############################################################################

if __name__ == '__main__':
    main()

##############################################################################

