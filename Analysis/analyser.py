#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 13:49:42 2017

@author: duman
"""

""" general simulation analyse module"""

##############################################################################

import sys
sys.path.append('../Utility')

import read_write

##############################################################################  
     
class Analyser:
    
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
        
    
    def write_results(self, write_func):
        """ write the analysis results"""
        
        write_func(self.results, self.savebase, \
                              self.analysisname, self.simdata)
        
        return
        
##############################################################################
