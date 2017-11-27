#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 11:10:46 2017

@author: duman
"""

""" non-gaussian parameter analysis"""

##############################################################################

import sys
sys.path.append('../include')

import argparse
import numpy as np

import read_write
import analyser

##############################################################################

def calculate_drift_subtracted_nongauss(polymers, sim):
    """ calculate and average the drift subtracted
    non-gaussian parameter (time and ensemble avgs)"""

    xu = polymers.xu
    ndelay = int(sim.nsteps/2)
    delay = np.zeros((ndelay), dtype=np.int64)
    alpha = np.zeros((ndelay), dtype=np.float64)

    for d in range(1, ndelay):
        delay[d] = d*sim.dt
        xd = xu[d:,0,:]-xu[:-d,0,:]
        xd_mean = np.mean(xd, axis=1)
        yd = xu[d:,1,:]-xu[:-d,1,:]
        yd_mean = np.mean(yd, axis=1)
        xd -= xd_mean[:,None]
        yd -= yd_mean[:, None]
        m2d = np.mean(np.mean( \
            xd**2 + yd**2, axis=1) )
        m4d = np.mean(np.mean( \
            xd**4 + yd**4, axis=1) )
        alpha[d] = m4d/(2*m2d**2) - 1

    return delay, alpha

##############################################################################

class AnalyseNonGaussian(analyser.Analyser):
    """ non-gaussian parameter analysis"""

    def __init__(self):

        args = self.get_args()
        analyser.Analyser.__init__(self, args.datafolder, args.savebase, \
                                   args.analysisname)
        self.read_data(args.beadsorpols, args.simtype)
        self.perform_analysis(calculate_drift_subtracted_nongauss)
        self.write_results(read_write.write_2d_analysis_data)

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

    nongauss_analysis = AnalyseNonGaussian()

    return

##############################################################################

if __name__ == "__main__":
    main()

##############################################################################

