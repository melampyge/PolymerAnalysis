#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 11:10:46 2017

@author: duman
"""

""" 4 point susceptibility analysis"""

##############################################################################

import sys
sys.path.append('../include')

import argparse
import numpy as np
import numpy.ma as ma

import read_write
import analyser

##############################################################################

def calc_4_pt_suscp(polymers, sim):
    """ calculate and average the susceptibility"""

    xu = polymers.xu
    threshold_amplitude = sim.avg_radius*2

    ndelay = int(sim.nsteps/2)
    delay = np.zeros((ndelay), dtype=np.int64)
    overlap = np.zeros((ndelay), dtype=np.float64)
    overlap_sq = np.zeros((ndelay), dtype=np.float64)
    suscp = np.zeros((ndelay), dtype=np.float64)
    overlap[0] = 1.0

    for d in range(1, ndelay):
        delay[d] = d*sim.dt

        ### subtract the mean displacement from cell displacements at the given delay time

        xd = xu[d:,0,:]-xu[:-d,0,:]
        xd_mean = np.mean(xd, axis=1)
        yd = xu[d:,1,:]-xu[:-d,1,:]
        yd_mean = np.mean(yd, axis=1)
        xd -= xd_mean[:,None]
        yd -= yd_mean[:,None]

        ### calculate the displacement magnitudes per cell per time

        displacements = np.sqrt(xd**2 + yd**2)

        for t0 in range(int(sim.nsteps-ndelay)):

            ### calculate the ensemble averaged overlap function

            masked = ma.masked_greater(displacements[t0], threshold_amplitude)  # filter large displacements
            masked = masked[~masked.mask]                                       # get the low displacements
            low_displacements = float(len(masked.data))                         # get the number of low displacements
            overlap_ens_avg = low_displacements / sim.npols                     # ens. averaged overlap fnc.
            overlap[d] += overlap_ens_avg
            overlap_sq[d] += overlap_ens_avg**2

        overlap[d] /= (sim.nsteps-ndelay)
        overlap_sq[d] /= (sim.nsteps-ndelay)
        suscp[d] = sim.npols*(overlap_sq[d] - overlap[d]**2)

    return delay, suscp

##############################################################################

class Analyse4ptsusc(analyser.Analyser):
    """ mean square displacement analysis"""

    def __init__(self):

        args = self.get_args()
        analyser.Analyser.__init__(self, args.datafolder, args.savebase, \
                                   args.analysisname)
        self.read_data(args.beadsorpols, args.simtype)
        self.perform_analysis(calc_4_pt_suscp)
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
                            help="Specific folder for saving, as in 4_pt_suscp")
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

    ssc_analysis = Analyse4ptsusc()

    return

##############################################################################

if __name__ == "__main__":
    main()

##############################################################################

