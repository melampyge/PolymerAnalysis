#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 11:10:46 2017

@author: duman
"""

""" mean square displacement analysis"""

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np

import read_write
import analyser
 
##############################################################################
    
def calculate_MSD(polymers, sim):
    """ calculate and average the mean squared displacement
    (time and ensemble avgs)"""
        
    xu = polymers.xu
    ndelay = int(sim.nsteps/2)
    delay = np.zeros((ndelay), dtype=np.int64)
    msd = np.zeros((ndelay), dtype=np.float64)
    
    for d in range(1, ndelay):
        delay[d] = d*sim.dt
        msd[d] = np.mean(np.mean( \
            (xu[d:,0,:]-xu[:-d,0,:])**2 + \
                (xu[d:,1,:]-xu[:-d,1,:])**2, axis=1) )
        
    return delay, msd      

##############################################################################
    
def calculate_drift_subtracted_MSD(polymers, sim):
    """ calculate and average the drift subtracted MSD (time and ensemble avgs)"""
        
    xu = polymers.xu
    ndelay = int(sim.nsteps/2)
    delay = np.zeros((ndelay), dtype=np.int64)
    msd = np.zeros((ndelay), dtype=np.float64)
    
    for d in range(1, ndelay):
        delay[d] = d*sim.dt
        xd = xu[d:,0,:]-xu[:-d,0,:]
        xd_mean = np.mean(xd, axis=1)
        yd = xu[d:,1,:]-xu[:-d,1,:]
        yd_mean = np.mean(yd, axis=1)
        xd -= xd_mean[:,None]
        yd -= yd_mean[:, None]
        msd[d] = np.mean(np.mean( \
            xd**2 + yd**2, axis=1) )
        
    return delay, msd  
    
##############################################################################
    
class AnalyseMSD(analyser.Analyser):
    """ mean square displacement analysis"""
    
    def __init__(self):
        
        args = self.get_args()
        analyser.Analyser.__init__(self, args.datafolder, args.savebase, \
                                   args.analysisname, args.choice)
        self.read_data(self.choice)
        self.perform_analysis(calculate_MSD)
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
        parser.add_argument("-c", "--choice", nargs="?", \
                            const = "simulation", \
                            help="Read choice --simulation, beads, polymers, or all--")        
        args = parser.parse_args()
        
        return args     
        
############################################################################## 
    