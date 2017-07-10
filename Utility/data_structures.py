
""" Data structures for storing information"""

##############################################################################

import numpy as np

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
    
    def __init__(self, datafolder, dt, density, nsteps, nbeads, \
                 npols, nbpp, bl, sigma, lx, ly, kappa, fp):
        
        Simulation.__init__(self, datafolder, dt, density, nsteps, nbeads, \
                 npols, nbpp, bl, sigma, lx, ly)
        self.kappa = kappa
        self.fp = fp
        
        ### define general physical parameters
        # NOTE THAT : 
        # these are defined from the longest length scale polymers
        
        self.length = self.get_length_of_polymer(self.nbpp[0])
        self.pe = self.get_pe_of_polymer(self.nbpp[0])
        self.xil = self.get_xil_of_polymer(self.nbpp[0])
        self.tau_diff = self.length**3 * self.gamma_0 * self.nbpp / self.length \
            / self.kT / 4.
        self.tau_advec = self.length * self.gamma_0 / self.fp
        self.vc = self.fp / self.gamma_0
        
        return
    
    def get_length_of_polymer(self, nb):
        """ calculate the length of the polymer"""
        
        return nb*self.bond_length
    
    def get_pe_of_polymer(self, nb):
        """ calculate the peclet number of the polymer"""
        
        l = self.get_length_of_polymer(nb)
        
        return self.fp*l**2/self
    
    def get_xil_of_polymer(self, nb):
        """ calculate the persistence length of the polymer"""
        
        l = self.get_length_of_polymer(nb)

        return self.kappa/self.kT/l
    
##############################################################################

class Subplots:
    """ plot structure"""
    
    totcnt = -1             # Total number of subplots 
    
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
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