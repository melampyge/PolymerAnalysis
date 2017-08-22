#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 17:05:32 2017

@author: duman
"""

""" Plot velocity structure factor"""

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

def five_thirds(x, a):

    return a * x**(-5./3.)

##############################################################################

def third(x, a):

    return a * x**(-3.)

##############################################################################

class PlotProps:
    """ encapsulate plot properties"""

    def __init__(self):

        self.read_fnc = read_write.read_2d_analysis_data
        self.num_ticks = 5
        self.ax_len = 1.0                          # Length of one subplot square box
        self.ax_b = 0.0                            # Beginning/offset of the subplot in the box
        self.ax_sep = 0.0                          # Separation length between two subplots
        self.total_subplots_in_x = 1               # Total number of subplots

        return

##############################################################################

def select_data(data, keys):
    """ select and/or normalize the data for plotting"""

    xd = {}
    yd = {}
    simd = {}
    for j, key in enumerate(keys):
        simd[key] = data[key].sim
        kxnorm = 2*np.pi/2/simd[key].avg_radius
        yd[key] = data[key].data[1]
        xd[key] = data[key].data[0]
        if type(xd[key]) != float:
            yd[key] /= xd[key]
        xd[key] /= kxnorm

    return xd, yd, simd

##############################################################################

@plotter.plot_2d
def plot_analysis(data, sep, savepdf, *args):
    """ plot the analysis data specific for this analysis"""

    fig = args[0]       # figure handle
    ax0 = args[1]       # axes handle

    ### normalize or select data

    keys = np.sort(data.keys())
    xd, yd, simd = select_data(data, keys)
    gsim = simd[simd.keys()[0]]
    avg_radius = gsim.avg_radius
    data_separator.assign_physicalvalues(sep.fixed_param, gsim)

    ### plot the data

    for key in keys:

        x = xd[key]
        y = yd[key]
        sim = simd[key]

        ### curve fitting

        if type(x) == float:
            continue
        popt, pcov = curve_fit(third, x[15:-10], y[15:-10])
        yfit = third(x[15:-10], popt[0])

        sep.legend_param.assign_physicalvalue(sim)
        label = data_separator.gen_label(sep.legend_param)
        line0 = ax0.plot(x, y, label=label, \
                         linewidth=2.0)
        line1 = ax0.plot(x[15:-10], yfit, '--', label='_nolegend_', \
                         linewidth=1.0, c='grey')

    ### scales

    ax0.set_xscale('log')
    ax0.set_yscale('log')

    ### title

    title = data_separator.gen_title(sep.fixed_param)
    ax0.set_title(title, fontsize=30)

    ### side labels

    ax0.set_xlabel(r'$k/k_{2R}$', fontsize=40)
    ax0.set_ylabel(r'$E(k)$', fontsize=40)

    ### limits

    #ax0.set_xlim((0, gsim.lx/2./avg_radius))
    #ax0.set_ylim(())

    ### ticks

    #ax0.xaxis.set_ticks(np.linspace(0, 15, num_ticks, endpoint=True))
    #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    ax0.tick_params(axis='both', which='major', labelsize=30)

    ### legend

    ax0.legend(bbox_to_anchor=(0.65, 0., 0.65, 1.), loc=2, borderaxespad=0., \
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

