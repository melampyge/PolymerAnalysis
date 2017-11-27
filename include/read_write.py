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
import data_separator

##############################################################################

def write_single_analysis_data(results, savebase, analysisname, sim):
    """ write single analysis data to the corresponding file"""

    ### create the path

    base = savebase + analysisname + '/'
    os.system("mkdir -p " + base)
    base += analysisname + '_'
    for j in range(len(sim.phase_params)):
        tag = data_separator.gen_tag(\
            sim.phase_params[j][0], sim.phase_params[j][1], "_")
        base += tag
        if j != len(sim.phase_params)-1:
            base += "_"
    fpath = base + '.txt'
    print fpath

    ### write the data

    fl = open(fpath, 'w')
    fl.write(str(results) + '\n')
    fl.close()

    return

##############################################################################

def write_1d_analysis_data(results, savebase, analysisname, sim):
    """ write 1d analysis data to the corresponding file"""

    ### create the path

    base = savebase + analysisname + '/'
    os.system("mkdir -p " + base)
    base += analysisname + '_'
    for j in range(len(sim.phase_params)):
        tag = data_separator.gen_tag(\
            sim.phase_params[j][0], sim.phase_params[j][1], "_")
        base += tag
        if j != len(sim.phase_params)-1:
            base += "_"
    fpath = base + '.txt'
    print fpath

    ### write the data

    fl = open(fpath, 'w')
    x = results
    N = len(x)
    for j in range(N):
        fl.write(str(x[j]) + '\n')

    fl.close()

    return

##############################################################################

def write_2d_analysis_data(results, savebase, analysisname, sim):
    """ write 2d analysis data to the corresponding file"""

    ### create the path

    base = savebase + analysisname + '/'
    os.system("mkdir -p " + base)
    base += analysisname + '_'
    for j in range(len(sim.phase_params)):
        tag = data_separator.gen_tag(\
            sim.phase_params[j][0], sim.phase_params[j][1], "_")
        base += tag
        if j != len(sim.phase_params)-1:
            base += "_"
    fpath = base + '.txt'
    print fpath

    ### write the data

    fl = open(fpath, 'w')
    x, y = results
    N = len(x)
    for j in range(N):
        fl.write(str(x[j]) + '\t\t' + str(y[j]) + '\n')

    fl.close()

    return

##############################################################################

def read_1d_analysis_data(f):
    """ read 1d analysis data"""

    if os.path.exists(f):
        data = np.transpose(np.loadtxt(f, dtype=np.float64))
        x = data
    else:
        x = 0.

    return x

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

