
""" A set of commonly used general helper functions"""

### 

##############################################################################

import numpy as np
import math

##############################################################################

def nearbyint(x):
    """ Round to the nearby integer"""
    
    if x >= 0:
        return math.floor(x+0.5)
    else:
        return math.floor(x-0.5)

##############################################################################
    
def min_img_dist(x1, x2, lx):
    """ compute the minimum image distance between two positions"""
    
    dx = x2 - x1 
    return dx-nearbyint(dx/lx)*lx

##############################################################################

def coords_per_pols(x, y, nbpp):
    """ get the coordinates of beads in terms of polymers"""
    
    splitter = np.cumsum(nbpp)[:-1]
    x_per_pol = np.split(x, splitter, axis=1)
    y_per_pol = np.split(y, splitter, axis=1)
    
    return x_per_pol, y_per_pol

##############################################################################

def gen_folder_path(base, d, k, f):
    """ generate folder path address from the parameters"""
    
    return base + 'density_' + str(d) + \
        '/kappa_' + str(k) + '/fp_' + str(f) + '/'

##############################################################################

def calc_velocities(x, d, dt):
    """ calculate the velocities"""
    
    return (x[d:, :, :] - x[:-d, :, :])/dt
    
##############################################################################

def get_img_pos(x, lx):
    """ get the image position in the central box 
    -- can be numpy array or single pos"""
    
    return x-np.floor(x/lx)*lx

##############################################################################

def get_pol_ids(nbeads, nbpp):
    """ get the polymer identities of beads"""
    
    pid = np.zeros((nbeads), dtype=np.int32)
    splitter = np.cumsum(nbpp)
    begin_idx = 0
    for j, end_idx in enumerate(splitter):
        print j, begin_idx, end_idx
        pid[begin_idx : end_idx] = j
        begin_idx = end_idx
    
    return
    
##############################################################################

def calc_bond_orientations(x, nsteps, nbeads, npols, nbpp):
    """ calculate the bond orientations"""
    
    ori = np.zeros((nsteps, nbeads), dtype=np.float32)
    for step in xrange(nsteps):
        k = 0
        for n in xrange(npols):
            for j in xrange(nbpp[n]-1):
                dr = x[step,:,k+1] - x[step,:,k]
                ori[step][k] = np.atan2(dr[1], dr[0])
                k += 1
            dr = x[step,:,k+1] - x[step,:,k]
            ori[step][k] = np.atan2(dr[1], dr[0])
            k += 1                
    
    return ori
    
##############################################################################

