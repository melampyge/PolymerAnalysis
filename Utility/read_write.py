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
    
def read_polymer_data(folder):
    """ read polymer data from hdf5 file"""
    
    ### file path
    
    fpath = folder + 'out.h5'
    assert os.path.exists(fpath), "The out.h5 file does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')
    
    ### polymer information
    
    comu = np.array(fl['/pols/comu'], dtype=np.float32)
        
    return data_structures.Polymers(comu)

##############################################################################
    
def read_bead_data(folder):
    """ read bead data from hdf5 file"""
    
    ### file path
    
    fpath = folder + 'out.h5'
    assert os.path.exists(fpath), "The out.h5 file does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')
    
    ### polymer information
    
    xu = np.array(fl['/beads/xu'], dtype=np.float32)
        
    return data_structures.Beads(xu)
    
##############################################################################
    
def read_sim_info(folder):
    """ read simulation info from hdf5 file,
    specific for bidisperse simulations"""
    
    ### file path
    
    fpath = folder + 'out.h5'
    assert os.path.exists(fpath), "The out.h5 file does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')    
    
    ### general information about the simulation
    
    nsteps = int(fl['/sim/nsteps'][...])
    nbeads = int(fl['/sim/nbeads'][...])
    npols = int(fl['/sim/npols'][...])
    nbpp = np.array(fl['/sim/nbpp'], dtype=np.int32)
    
    ### simulation parameters
    
    density = float(fl['/params/density'][...])
    kappa = float(fl['/params/kappa'][...])
    fp = float(fl['/params/fp'][...])
#    bl = float(fl['/params/bl'][...])
#    sigma = float(fl['/params/sigma'][...])
#    dt = float(fl['/params/dt'][...])
#    lx = float(fl['/params/lx'][...])
#    ly = float(fl['/params/ly'][...])
    bl = 0.5
    sigma = 1.0
    dt = 50.0
    lx = 221.6
    ly = 221.6
    
    fl.close()
    
    ### generate classes to submerge data
    
    sim = data_structures.SimulationBidispersePolymers(folder, \
                                     dt, density, \
                                     nsteps, nbeads, \
                                     npols, nbpp, bl, sigma, \
                                     lx, ly, kappa, fp)
    
    return sim
    
##############################################################################

def read_h5_file(folder):
    """ read all of the hdf5 file"""
    
    beads = read_bead_data(folder)
    pols = read_polymer_data(folder)
    sim = read_sim_info(folder)
    
    return beads, pols, sim
    
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
    
    data = np.transpose(np.loadtxt(f, dtype=np.float64))
    x = data[0]
    y = data[1]

    return x, y   
    
##############################################################################
    