#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 17:05:32 2017

@author: duman
"""

""" Plot mean square displacement function"""

###

##############################################################################

import sys
sys.path.append('../../include')

import os
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

def lin(x, a):

    return a * x

##############################################################################

def quad(x, a):

    return a * x**2

##############################################################################

def th_msd(t, Dt, Dr, v0):

    return (4*Dt+v0**2/Dr)*t + v0**2*(np.exp(-2*t*Dr-1))/(2*Dr**2)

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
        yd[key] = data[key].data[1]/simd[key].length**2
        xd[key] = data[key].data[0]/simd[key].tau_diff
        if type(data[key].data[0]) == float:
            xd[key] = np.ones(30)*xd[key]
            yd[key] = np.ones(30)*yd[key]

    return xd, yd, simd

##############################################################################

def calc_exponent(x, y):
    """ calculate the exponent from MSD data"""

    xcut = x[70:]
    ycut = y[70:]

    #s = interpolate.UnivariateSpline(xcut, ycut)
    #xnew = np.linspace(np.min(xcut), np.max(xcut), num=100*len(xcut))
    #xnew = np.logspace(np.min(xcut), np.max(xcut), num=10*len(xcut), base=10.0)
    #ynew = s(xnew)
    #intp = interpolate.interp1d(xnew, ynew)
    #ynew = intp(xnew)

    xnew = np.log10(xcut)
    ynew = np.log10(ycut)

    return (ynew[-1]-ynew[0])/(xnew[-1]-xnew[0])

##############################################################################

def extract_exponent(xd, yd, simd, keys):
    """ extract and write the exponent"""

    base = "/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/Exponents/"
    os.system("mkdir -p " + base)
    base += "Exponents_"
    for j, key in enumerate(keys):
        sim = simd[key]
        x = xd[key]
        y = yd[key]

        if sim.pe is -1:
            continue

        val = calc_exponent(x, y)
        print val
        path = "eps_" + str(sim.eps) + "_fp_" + str(sim.fp) + \
            "_areak_" + str(sim.areak) + "_kappa_" + str(sim.kappa) + ".txt"
        flpath = base + path
        print flpath

        fl = open(flpath, "w")
        fl.write(str(val) + "\n")
        fl.close()

    return

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
    avg_radius = gsim.length
    data_separator.assign_physicalvalues(sep.fixed_param, gsim)

    ### perform extra analysis

    #extract_exponent(xd, yd, simd, keys)

    ### plot the data

    cnt = 0
    for key in keys:

        x = xd[key][2:]
        y = yd[key][2:]
        sim = simd[key]

        ### curve fitting

        popt, pcov = curve_fit(lin, x, y)
        yfit_lin = lin(x, popt[0])

        popt, pcov = curve_fit(quad, x, y)
        yfit_quad = quad(x, popt[0])

        #yth = th_msd(x*sim.tau_diff, 1.0/40.0, 0.005, sim.fp)

        sep.legend_param.assign_physicalvalue(sim)
        label = data_separator.gen_label(sep.legend_param)
        line0 = ax0.plot(x, y, label=label, \
                         linewidth=3.0)

        if cnt == 0:
            line1 = ax0.plot(x, yfit_lin, "--", c='grey', \
                             label="_nolegend_", linewidth=2.0)
            line2 = ax0.plot(x, yfit_quad, "--", c='grey', \
                             label="_nolegend_", linewidth=2.0)
            #line3 = ax0.plot(x, yth, "--", c='grey', \
            #                 label="_nolegend_", linewidth=2.0)
            cnt += 1

    ### scales

    ax0.set_xscale('log')
    ax0.set_yscale('log')

    ### title

    title = data_separator.gen_title(sep.fixed_param)
    ax0.set_title(title, fontsize=20)

    ### side labels

    ax0.set_xlabel(r'$t/\tau_{D}$', fontsize=30)
    ax0.set_ylabel(r'$\Delta r^2/L^2$', fontsize=30)

    ### limits

    #ax0.set_xlim((1, 700.0))
    #ax0.set_xlim((0, gsim.lx/2./avg_radius))
    #ax0.set_ylim(())

    ### ticks

    #ax0.xaxis.set_ticks(np.linspace(0, 15, num_ticks, endpoint=True))
    #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    ax0.tick_params(axis='both', which='major', labelsize=20)

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

