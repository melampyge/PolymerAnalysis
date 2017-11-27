#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 11:10:46 2017

@author: duman
"""

""" overlap function analysis"""

##############################################################################

import sys
sys.path.append('../include')

import argparse
import numpy as np
import numpy.ma as ma

import read_write
import analyser

##############################################################################

def calculate_overlap_fnc(polymers, sim):
    """ calculate and average the overlap fnc (time and ensemble avgs)"""

    xu = polymers.xu
    threshold_amplitude = sim.avg_radius*2

    ndelay = int(sim.nsteps/2)
    delay = np.zeros((ndelay), dtype=np.int64)
    overlap = np.zeros((ndelay), dtype=np.float64)
    overlap[0] = 1.0

    for d in range(1, ndelay):
        delay[d] = d*sim.dt

        ### subtract the mean displacement from cell displacements at the given delay time

        xd = xu[d:,0,:]-xu[:-d,0,:]
        xd_mean = np.mean(xd, axis=1)
        yd = xu[d:,1,:]-xu[:-d,1,:]
        yd_mean = np.mean(yd, axis=1)
        xd -= xd_mean[:,None]
        yd -= yd_mean[:, None]

        ### calculate the displacement magnitudes per particle

        displacements = np.sqrt(xd**2 + yd**2)

        ### calculate the overlap function by taking an ensemble and time average

        masked = ma.masked_greater(displacements, threshold_amplitude)
        masked = masked[~masked.mask]
        low_displacements = float(len(masked.data))
        total_displacements = float(np.size(displacements))
        overlap[d] = low_displacements / total_displacements

    return delay, overlap

##############################################################################

class AnalyseOverlapFnc(analyser.Analyser):
    """ mean square displacement analysis"""

    def __init__(self):

        args = self.get_args()
        analyser.Analyser.__init__(self, args.datafolder, args.savebase, \
                                   args.analysisname)
        self.read_data(args.beadsorpols, args.simtype)
        self.perform_analysis(calculate_overlap_fnc)
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
                            help="Specific folder for saving, as in Overlap_fnc")
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

    overlap_analysis = AnalyseOverlapFnc()

    return

##############################################################################

if __name__ == "__main__":
    main()

##############################################################################

