#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 16:31:55 2017

@author: duman
"""

""" Data types and functions for plotting"""

### 

##############################################################################

import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import read_write
import misc_tools 

import seaborn as sns
sns.set(style="white",context='paper',
        font_scale=1.2,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 30, 
            'font.family': 'sans',"figure.dpi":300,
            "xtick.major.size": 8, "ytick.major.size": 8,
            'grid.linestyle': '--'})   

##############################################################################

class Subplots:
    """ subplots structure"""
    
    totcnt = -1             # Total number of subplots -- static member
    
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
        return
        
    def addSubplot(self):
        """ add a subplot in the grid structure"""
        
        ## increase the number of subplots in the figure
        
        self.totcnt += 1
        
        ## get indices of the subplot in the figure
        
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])

##############################################################################

class GeneralPlot(Subplots):
    """ general plot structure"""
    
    def __init__(self, num_ticks=5, 
                 ax_len=1.0, ax_b=0.0, 
                 ax_sep=0.0, total_subplots_in_x=1):
    
        self.fig = plt.figure()
        Subplots.__init__(self, self.fig, \
                          ax_len, ax_sep, ax_b, total_subplots_in_x)
        self.ax0 = self.addSubplot()
        
        return
        
##############################################################################    