#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 11:10:46 2017

@author: duman
"""

""" spiral number analysis"""

##############################################################################

import sys
sys.path.append('../include')

import argparse
import numpy as np
import math

import read_write
import analyser

##############################################################################

def calc_spiral_num_per_fil(x, y, nbpp):
    """ calculate the spiral number per filament"""

    phi = np.zeros((nbpp-1))
    for j in range(1, nbpp):
        dx = x[j]-x[j-1]
        dy = y[j]-y[j-1]
        dphi = math.atan2(dy, dx)
        phi[j-1] = dphi

    phi2 = np.copy(phi)
    nbonds = len(phi)
    for j in range(1, nbonds):
        dphi = phi[j] - phi[j-1]
        if dphi < -np.pi:
            dphi += 2*np.pi
        elif dphi > np.pi:
            dphi -= 2*np.pi
        phi2[j] = phi2[j-1] + dphi

    return (phi2[-1] - phi2[0])/2/np.pi

##############################################################################

def calculate_spiral_number(beads, sim):
    """ calculate and average the spiral number
    (time and ensemble avgs)"""

    spiral_number = np.zeros((sim.nsteps), dtype=np.float64)
    xu = beads.xu

    for step in range(sim.nsteps):
        k = 0
        s_per_step = 0.0
        for n in range(sim.npols):
            s_per_step += calc_spiral_num_per_fil(xu[step,0,k:k+sim.nbpp[n]],\
                xu[step,1,k:k+sim.nbpp[n]], sim.nbpp[n])
            k += sim.nbpp[n]
        spiral_number[step] = s_per_step/sim.npols

    glob_spiral_number = np.mean(np.abs(spiral_number))

    return glob_spiral_number

##############################################################################

class AnalyseSpiralNumber(analyser.Analyser):
    """ spiral number analysis"""

    def __init__(self):

        args = self.get_args()
        analyser.Analyser.__init__(self, args.datafolder, args.savebase, \
                                   args.analysisname)
        self.read_data(args.beadsorpols, args.simtype)
        self.perform_analysis(calculate_spiral_number)
        self.write_results(read_write.write_single_analysis_data)

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
                            help="Specific folder for saving, as in Spiral_number")
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

    spiral_analysis = AnalyseSpiralNumber()

    return

##############################################################################

if __name__ == "__main__":
    main()

##############################################################################

