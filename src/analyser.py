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
import data_structures

##############################################################################

class Analyser:

    def __init__(self, fd, sb, an):

        self.datafolder = fd            # folder containing simulation data
                                        # (h5 file)
        self.savebase = sb              # folder inside which analysis data
                                        # will be saved
        self.analysisname = an          # name of the analysis folder and file

        return


    def read_data(self, bop, sim_type):
        """ read simulation, polymers, and/or beads data"""

        self.read_choice = bop
        self.sim_type = sim_type

        if self.read_choice == "beads":
            self.beaddata = data_structures.Beads(self.datafolder)
        elif self.read_choice == "polymers":
            self.poldata = data_structures.Polymers(self.datafolder, self.sim_type)
        else:
            raise "The beadsorpols should be either polymers, or beads"

        if self.sim_type == "filaments":
            self.simdata = data_structures.SimulationFilaments(self.datafolder)
        elif self.sim_type == "cells":
            self.simdata = data_structures.SimulationCells(self.datafolder)
        else:
            raise "The simulation type should be either filaments, or cells"

        return


    def perform_analysis(self, analysis_func):
        """ perform the analysis on the selected data"""

        if self.read_choice == "beads":
            self.results = analysis_func(self.beaddata, self.simdata)
        elif self.read_choice == "polymers":
            self.results = analysis_func(self.poldata, self.simdata)
        else:
            self.results = analysis_func(self.simdata)

        return


    def write_results(self, write_func):
        """ write the analysis results"""

        write_func(self.results, self.savebase, \
                              self.analysisname, self.simdata)

        return

##############################################################################
