#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:15:10 2017

@author: duman
"""

""" data read and write functionality"""

##############################################################################

import numpy as np
import os
import h5py
import misc_tools
import data_structures

##############################################################################

def write_2d_analysis_data(results, savebase, analysisname, sim):
    """ write 2d analysis data to the corresponding file"""

    ### create the path

    base = savebase + analysisname + '/'
    os.system("mkdir -p " + base)
    base += analysisname + '_'
    folderpath = misc_tools.gen_folder_path(base, '_', sim.phase_params)
    fpath = folderpath + '.txt'

    ### write the data

    fl = open(fpath, 'w')
    x, y = results
    N = len(x)
    for j in range(N):
        fl.write(str(x[j]) + '\t\t' + str(y[j]) + '\n')

    fl.close()

    return

##############################################################################

def read_2d_analysis_data(f):
    """ read 2d analysis data"""

    if os.path.exists(f):
        data = np.transpose(np.loadtxt(f, dtype=np.float64))
        x = data[0]
        y = data[1]
    else:
        x = 0.
        y = 0.

    return x, y

##############################################################################

def read_single_analysis_data(f):
    """ read single analysis data"""

    if os.path.exists(f):
        data = np.loadtxt(f, dtype=np.float64)
    else:
        data = 0.

    return data

##############################################################################

def read_8d_analysis_data(f):
    """ read 8d analysis data"""

    if os.path.exists(f):
        data = np.transpose(np.loadtxt(f, dtype=np.float64))
    else:
        data = 0.

    return data

##############################################################################

