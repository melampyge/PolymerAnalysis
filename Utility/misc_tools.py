
""" A set of commonly used general helper functions"""

### 

##############################################################################

import numpy as np

##############################################################################

def gen_folder_path(base, d, k, f):
    """ generate folder path address from the parameters"""
    
    return base + 'density_' + str(d) + \
        '/kappa_' + str(k) + '/fp_' + str(f) + '/'

##############################################################################


