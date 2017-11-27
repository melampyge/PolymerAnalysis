#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 11:10:46 2017

@author: duman
"""

""" displacement distribtuion analysis"""

##############################################################################

import sys
sys.path.append('../include')

import argparse
import numpy as np

import read_write
import analyser

##############################################################################

def calculate_drift_subtracted_displacements(polymers, sim):
    """ calculate the drift subtracted displacements"""

    xu = polymers.xu
    d = 5

    xd = xu[d:,0,:]-xu[:-d,0,:]
    xd_mean = np.mean(xd, axis=1)
    yd = xu[d:,1,:]-xu[:-d,1,:]
    yd_mean = np.mean(yd, axis=1)
    xd -= xd_mean[:,None]
    yd -= yd_mean[:, None]
    disp = xd**2 + yd**2
    disp = disp.flatten()

    return disp

##############################################################################

class AnalyseDisplacement(analyser.Analyser):
    """ displacement distribution analysis"""

    def __init__(self):

        args = self.get_args()
        analyser.Analyser.__init__(self, args.datafolder, args.savebase, \
                                   args.analysisname)
        self.read_data(args.beadsorpols, args.simtype)
        self.perform_analysis(calculate_drift_subtracted_displacements)
        self.write_results(read_write.write_1d_analysis_data)

        return


    def get_args(self):
        """ get the command line arguments"""

        parser = argparse.ArgumentParser()
        parser.add_argument("-fd", "--datafolder", \
                            help="Folder containing data, as in /local/duman/SIMULATIONS/Bidisperse_Filaments/.../")
        parser.add_argument("-sb", "--savebase", nargs="?", \
                            const = "/usr/users/iff_th2/duman/Bidisperse_Filaments/DATA/", \
                            help="Folder to save the data, as in /usr/users/iff_th2/duman/Bidisperse_Filaments/DATA/")
        parser.add_argument("-an", "--analysisname", \
                            help="Specific folder for saving, as in MSD")
        parser.add_argument("-bop", "--beadsorpols", nargs="?", \
                            const = "beads", \
                            help="Read choice --beads or polymers--")
        parser.add_argument("-st", "--simtype", nargs="?", \
                            const = "simulation", \
                            help="Choose simulation type --filaments or cells--")
        args = parser.parse_args()

        return args

##############################################################################

def main():

    displacement_analysis = AnalyseDisplacement()

    return

##############################################################################

if __name__ == "__main__":
    main()

##############################################################################

