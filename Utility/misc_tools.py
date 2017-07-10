
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
    x_per_pol = np.split(x, splitter)
    y_per_pol = np.split(y, splitter)
    
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

def calc_bond_orientations(x, y, nbpp):
    """ calculate the bond orientations"""
    
    x_per_pol, y_per_pol = coords_per_pols(x, y, nbpp)
    dx = [xpp[1:]-xpp[:-1] for xpp in x_per_pol]
    dy = [ypp[1:]-ypp[:-1] for ypp in y_per_pol]
    print np.shape(dx)
    print dx
    exit(1)
    
    return
    
##############################################################################

