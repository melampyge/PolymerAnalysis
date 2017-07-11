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
    
    nsteps = fl['/sim/nsteps'][...]
    nbeads = fl['/sim/nbeads'][...]
    npols = fl['/sim/npols'][...]
    nbpp = np.array(fl['/sim/nbpp'], dtype=np.float32)
    
    ### simulation parameters
    
    density = fl['/params/density'][...]
    kappa = fl['/params/kappa'][...]
    fp = fl['/params/fp'][...]
#    bl = fl['/params/bl'][...]
#    sigma = fl['/params/sigma'][...]
#    dt = fl['/params/dt'][...]
#    lx = fl['/params/lx'][...]
#    ly = fl['/params/ly'][...]
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