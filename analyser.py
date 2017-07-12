#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 13:49:42 2017

@author: duman
"""

""" general simulation analyse module"""

##############################################################################

import sys
sys.path.append('Utility')
sys.path.append('Analysis')

import argparse
import numpy as np
import read_write

##############################################################################  

def get_args():
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
    parser.add_argument("-af", "--analysisfunction", \
                        help="Name of the analysis function, as in calculate_MSD")  
    parser.add_argument("-wf", "--writefunction", \
                        help="Name of the write function, as in write_2d_data")        
    args = parser.parse_args()
    
    return args

##############################################################################  
     
class Analyse:
    
    def __init__(self, fd, sb, an):
        
        self.datafolder = fd            # folder containing simulation data 
                                        # (h5 file)
        self.savebase = sb              # folder inside which analysis data 
                                        # will be saved
        self.analysisname = an          # name of the analysis folder and file
        
        return
        
        
    def read_data(self, choice):
        """ read simulation, polymers, and/or beads data"""
        
        self.read_choice = choice
        
        if self.read_choice == "beads":
            self.beaddata = read_write.read_bead_data(self.datafolder)
            self.simdata = read_write.read_sim_info(self.datafolder)
        elif self.read_choice == "polymers":
            self.poldata = read_write.read_polymer_data(self.datafolder)
            self.simdata = read_write.read_sim_info(self.datafolder)            
        elif self.read_choice == "sim":
            self.simdata = read_write.read_sim_info(self.datafolder)
        elif self.read_choice == "all":
            self.beaddata, self.poldata, self.simdata = \
                read_write.read_h5_file(self.datafolder)
        else:
            raise "The read choices should be either polymers, beads, simulation or all"
    
        return 
        
        
    def perform_analysis(self, analysis_func):
        """ perform the analysis on the selected data"""
        
        if self.read_choice == "beads":
            self.results = analysis_func(self.beaddata, self.simdata)
        elif self.read_choice == "polymers":
            self.results = analysis_func(self.poldata, self.simdata)
        elif self.read_choice == "sim":
            self.results = analysis_func(self.simdata)
        elif self.read_choice == "all":
            self.results = analysis_func(self.beaddata, self.poldata, \
                                         self.simdata)
        else:
            raise "The read choices should be either polymers, beads, simulation or all"    
            
        return
        
    
    def write_results(self, func):
        """ write the analysis results"""
        
        read_write.func(self.results, self.savebase, \
                        self.analysisname, self.simdata)
        
        return
        
##############################################################################

def main():
    
    ### get the command line arguments
    
    args = get_args()
    
    ### create an instance of specific data analysis 
    
    analyser = Analyse(args.datafolder, args.savebase, args.analysisname)
    
    ### read simulation data
    
    analyser.read_data(args.choice)
    
    ### perform the selected analysis
    
    analyser.perform_analysis(args.analysis_function)
    
    ### write the analysis results
    
    analyser.write_results(args.writefunction)
    
    return
   
##############################################################################
    
if __name__ == "__main__":
    main()

##############################################################################