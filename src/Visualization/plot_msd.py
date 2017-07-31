#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 17:05:32 2017

@author: duman
"""

""" Plot mean square displacement"""

###

##############################################################################

import sys
sys.path.append('../../include')

import argparse
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt

import plotter
import read_write
import misc_tools
import data_separator

##############################################################################

def get_args():
    """ get the command line arguments"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--density", type=float, nargs="?", const=-1, \
                        help="Packing fraction")
    parser.add_argument("-k", "--kappa", type=float, nargs="?", const=-1, \
                        help="Bending rigidity")
    parser.add_argument("-f", "--fp", type=float, nargs="?", const=-1, \
                        help="Propulsion force")
    parser.add_argument("-e", "--eps", type=float, nargs="?", const=-1, \
                        help="Strength of LJ potential")
    parser.add_argument("-a", "--areak", type=float, nargs="?", const=-1, \
                        help="Area constraint potential strength")
    parser.add_argument("-fd", "--folder", \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Bidisperse_Filaments/Simulations/")
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Bidisperse_Filaments/PLOTS/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Bidisperse_Filaments/PLOTS/")
    parser.add_argument("-ab", "--analysisbase", nargs="?", \
                        const="/usr/users/iff_th2/duman/Bidisperse_Filaments/DATA/", \
                        help="Folder to load the data from, as in /usr/users/iff_th2/duman/Bidisperse_Filaments/DATA/")
    parser.add_argument("-an", "--analysisname", nargs="?", \
                        const="MSD", \
                        help="Name of the analysis, as in MSD")
    parser.add_argument("-s","--savepdf", action="store_true", help="Decide whether to save in pdf or not")
    args = parser.parse_args()

    if type(args.eps) and type(args.areak) == None:
        legend_param, fixed_params = data_separator.read_highdens_fil_param_data(args)
        args.simtype = "filaments"
    else:
        legend_param, fixed_params = data_separator.read_cell_param_data(args)
        args.simtype = "cells"

    print args

    return args, legend_param, fixed_params

##############################################################################

def set_plot_props(args, param_choice, fixed_params):
    """ set plot properties"""

    xlab = r"\Delta t/tau_D"
    ylab = r"\Delta_r^2/L^2"
    title = ""
    legend = ""

    return xlab, ylab, title, legend

##############################################################################

def normalize_data(x, data, sims):
    """ normalize the axis"""

    for j, key in enumerate(sims.keys()):
        x[key] = x[key]/sims[key].tau_diff
        data[key] = data[key]/sims[key].length**2

    return

##############################################################################

def main():

    args, legend_param, fixed_params = get_args()
    x, data, param_choice, sims = data_separator.rearrange_cell_data(args, \
        legend_param, fixed_params, read_write.read_2d_analysis_data)
    p = plotter.GeneralPlot()
    xlab, ylab, title, legend = set_plot_props(args, param_choice, fixed_params)
    normalize_data(x, data, sims)
    p.plot_2d(x, data, sims, args.savebase, args.analysisname, \
              xlab, ylab, title, legend, param_choice, args.savepdf)

    return

##############################################################################

if __name__ == '__main__':
    main()

##############################################################################
