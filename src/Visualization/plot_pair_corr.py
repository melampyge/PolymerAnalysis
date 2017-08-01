#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 17:05:32 2017

@author: duman
"""

""" Plot pair correlation function"""

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

def get_args(read_fnc):
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
                        const="Pair_corr", \
                        help="Name of the analysis, as in Pair_corr")
    parser.add_argument("-st", "--simtype", \
                        help="Simulation type --should be highdensfilaments, bidispersefilaments or cells--")
    parser.add_argument("-s","--savepdf", action="store_true", help="Decide whether to save in pdf or not")
    args = parser.parse_args()

    separator = data_separator.Separator(args, read_fnc)

    return separator, args.folder, args.savebase, args.analysisbase, \
        args.analysisname, args.savepdf

##############################################################################

class PlotProps:
    """ encapsulate plot properties"""

    def __init__(self):

        self.xlab = r"$\Delta r/2R$"
        self.ylab = r"$g(r)$"
        self.title = ""
        self.legend = ""
        self.xscale = "linear"
        self.yscale = "linear"
        self.set_xlim = True
        self.set_ylim = False
        self.xlim = (0., 10.)
        self.ylim = (0., 100.)
        self.set_xticks = False
        self.set_yticks = False
        self.xticks = []
        self.yticks = []

        return

##############################################################################

def normalize_data(x, data, sims):
    """ normalize the axis"""

    for j, key in enumerate(sims.keys()):
        x[key] = x[key]/2./sims[key].avg_radius

    return

##############################################################################

def main():

    separator, folder, savebase, analysisbase, analysisname, savepdf = get_args(\
        read_write.read_2d_analysis_data)
    x, data, param_choice, sims = separator.out
    p = plotter.GeneralPlot()
    plotprops = PlotProps()
    normalize_data(x, data, sims)
    p.plot_2d(x, data, sims, savebase, analysisname, \
              param_choice, separator, \
              plotprops, savepdf)

    return

##############################################################################

if __name__ == '__main__':
    main()

##############################################################################
