
""" Data structures for storing information"""

##############################################################################

import numpy as np
from collections import OrderedDict
import os 
import h5py

##############################################################################
        
class Beads:
    """ container for bead information"""
    
    def __init__(self, folder):
        
        self.read_data(folder)
        
        return

    def read_data(self, folder):
        """ read bead data from hdf5 file"""
        
        ### file path
        
        fpath = folder + 'out.h5'
        assert os.path.exists(fpath), "The out.h5 file does NOT exist for "\
             + fpath
        fl = h5py.File(fpath, 'r')
        
        ### bead information
        
        self.xu = np.array(fl['/pols/comu'], dtype=np.float32)
            
        return 
        
    def get_pol_ids(self, nbeads, nbpp):
        """ get the polymer identities of beads"""
        
        self.pid = np.zeros((nbeads), dtype=np.int32)
        splitter = np.cumsum(nbpp)
        begin_idx = 0
        for j, end_idx in enumerate(splitter):
            self.pid[int(begin_idx) : int(end_idx)] = j
            begin_idx = end_idx
        
        return   

        
    def calc_bond_orientations(self, x, nsteps, nbeads, npols, nbpp):
        """ calculate the bond orientations"""
        
        self.ori = np.zeros((nsteps, nbeads), dtype=np.float32)
        for step in xrange(nsteps):
            k = 0
            for n in xrange(npols):
                for j in xrange(nbpp[n]-1):
                    dr = x[step,:,k+1] - x[step,:,k]
                    self.ori[step][k] = np.atan2(dr[1], dr[0])
                    k += 1
                dr = x[step,:,k+1] - x[step,:,k]
                self.ori[step][k] = np.atan2(dr[1], dr[0])
                k += 1                
        
        return 
        
##############################################################################

class Polymers:
    """ container for polymer information"""
    
    def __init__(self, folder):
        
        self.read_data(folder)

        return
        
    def read_data(self, folder):
        """ read polymer data from hdf5 file"""
        
        ### file path
        
        fpath = folder + 'out.h5'
        assert os.path.exists(fpath), "The out.h5 file does NOT exist for "\
             + fpath
        fl = h5py.File(fpath, 'r')
        
        ### polymer information
        
        self.xu = np.array(fl['/pols/comu'], dtype=np.float32)
            
        return       
          
##############################################################################

class Simulation:
    """ container for general simulation information"""
    
    def __init__(self, datafolder, dt, density, nsteps, nbeads, \
                 npols, nbpp, bl, sigma, lx, ly):
        
        self.folder = datafolder        # data folder
        self.dt = dt                    # timestep between two data points
        self.density = density          # packing fraction
        self.nsteps = nsteps            # number of data points
        self.nbeads = nbeads            # number of beads
        self.npols = npols              # number of polymers
        self.nbpp = nbpp                # number of beads per polymer (descending sorted)
        self.bond_length = bl           # bond length
        self.sigma = sigma              # effective bead radius (assuming LJ)
        self.lx = lx                    # box length in x 
        self.ly = ly                    # box length in y  
        self.kT = 1.                    # kb * T
        self.gamma_0 = 1.               # friction coefficient per bead
            
        return

##############################################################################

class SimulationBidispersePolymers(Simulation):
    """ container for general simulation information for bidisperse polymers"""
    
    def __init__(self, datafolder):
        
        self.read_sim_info(datafolder)
        Simulation.__init__(self, datafolder, self.dt, \
                            self.density, self.nsteps, self.nbeads, \
                            self.npols, self.nbpp, self.bl, \
                            self.sigma, self.lx, self.ly)
        
        ### define general physical parameters
        # NOTE THAT : 
        # these are defined from the longest length scale polymers
        
        self.length = self.get_length_of_polymer(self.nbpp[0])
        self.pe = self.get_pe_of_polymer(self.nbpp[0])
        self.xil = self.get_xil_of_polymer(self.nbpp[0])
        self.tau_diff = self.length**3 * self.gamma_0 * self.nbpp[0] / self.length \
            / self.kT / 4.
        self.tau_advec = self.length * self.gamma_0 / self.fp
        self.vc = self.fp / self.gamma_0
        
        ### generate key parameters that map the phase space
        
        self.gen_key_params()
        
        return
    
    def get_length_of_polymer(self, nb):
        """ calculate the length of the polymer"""
        
        return nb*self.bond_length
    
    def get_pe_of_polymer(self, nb):
        """ calculate the peclet number of the polymer"""
        
        l = self.get_length_of_polymer(nb)
        
        return self.fp*l**2/self.kT
    
    def get_xil_of_polymer(self, nb):
        """ calculate the persistence length of the polymer"""
        
        l = self.get_length_of_polymer(nb)

        return self.kappa/self.kT/l
        
    def gen_key_params(self):
        """ generate key parameters that help the mapping of the phase space"""
        
        key_params = ([('density', self.density), \
                       ('kappa', self.kappa), ('fp', self.fp)])
        self.phase_params = key_params
        
        return
 
    def read_sim_info(self, folder):
        """ read simulation info from hdf5 file,
        specific for bidisperse simulations"""
        
        ### file path
        
        fpath = folder + 'out.h5'
        assert os.path.exists(fpath), "The out.h5 file does NOT exist for " + fpath
        fl = h5py.File(fpath, 'r')    
        
        ### general information about the simulation
        
        self.nsteps = int(fl['/sim/nsteps'][...])
        self.nbeads = int(fl['/sim/nbeads'][...])
        self.npols = int(fl['/sim/npols'][...])
        self.nbpp = np.array(fl['/sim/nbpp'], dtype=np.int32)
        
        ### simulation parameters
        
        self.density = float(fl['/params/density'][...])
        self.kappa = float(fl['/params/kappa'][...])
        self.fp = float(fl['/params/fp'][...])
    #    self.bl = float(fl['/params/bl'][...])
    #    self.sigma = float(fl['/params/sigma'][...])
    #    self.dt = float(fl['/params/dt'][...])
    #    self.lx = float(fl['/params/lx'][...])
    #    self.ly = float(fl['/params/ly'][...])
        self.bl = 0.5
        self.sigma = 1.0
        self.dt = 50.0
        self.lx = 221.6
        self.ly = 221.6
        
        fl.close()
        
        return
        
##############################################################################
     