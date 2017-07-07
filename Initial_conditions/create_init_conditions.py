# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 16:42:23 2017

@author: duman
"""

""" Create initial conditions 
    for a bidisperse ensemble of filaments"""
    
### NOTE THAT:
# I am adopting the notiation of polymers to generalize cells and filaments
# in terms of lingo
    
###########################################################################

import numpy as np
import argparse
import random
import math

###########################################################################

class Simulation:   
    
    
    def __init__(self, llong, lshort, dens, npl):
        
        ### set prefixed properties within the script
        
        self.set_prefixed_props()
        
        ### set command line input properties selected
        
        self.set_inputopt_props(llong, lshort, dens, npl)
        
        ### infer general polymer information 
        
        self.gen_pol_info()
        
        ### infer information about the simulation box
        
        self.gen_box_info()
        
        return
        
        
    def set_prefixed_props(self):
        
        self.bl = 0.5                   # bond length
        self.sigma = 1.0                # LJ radius
        
        return 
        
        
    def set_inputopt_props(self, llong, lshort, dens, npols_l):
        
        lpols = {}
        lpols['long'] = llong
        lpols['short'] = lshort
        self.lpols = lpols              # polymer lengths
        
        npols = {}
        npols['long'] = npols_l         
        self.npols = npols              # number of polymers
        
        self.density = dens             # packing fraction
        
        return
        
        
    def gen_pol_info(self):
        
        nbpp = {}
        nbpp['long'] = int(self.lpols['long']/self.bl+1)
        nbpp['short'] = int(self.lpols['short']/self.bl+1)
        self.nbeads_per_pol = nbpp      # number of beads per polymer
        
        nbeads = {}
        nbeads['long'] = self.npols['long']*self.nbeads_per_pol['long']
        nbeads['short'] = nbeads['long']
        self.nbeads = nbeads            # total number of beads per polymer type
        
        self.totnbeads = self.nbeads['long'] + self.nbeads['short']
                                        # total number of beads
                                        
        self.npols['short'] = int(self.nbeads['short']/self.nbeads_per_pol['short'])
                                
        
        self.totnpols = self.npols['long'] + self.npols['short']
        
        return
        
        
    def gen_box_info(self):
        
        area_of_pol = np.sum(\
            [self.lpols[k]*self.sigma*self.npols[k] for k in self.lpols.keys()])
        self.area_box = area_of_pol/self.density
        self.box_length = np.sqrt(self.area_box)
        
        return
        
        
    def gen_polymers(self):
        
        random.seed()   # NOTE: the scope of PRNG is all the file!
        
        pols = []
        pol_id = 0
        bead_id = 0
        for j in xrange(self.npols['long']):
            pols.append( Polymer(pol_id+1, self.nbeads_per_pol['long'], \
                self.bl, self.box_length, bead_id) )
            pol_id += 1
            bead_id += self.nbeads_per_pol['long']
            
        for j in xrange(self.npols['short']):
            pols.append( Polymer(pol_id+1, self.nbeads_per_pol['short'], \
                self.bl, self.box_length, bead_id) )
            pol_id += 1
            bead_id += self.nbeads_per_pol['short']  
          
        return pols
   
###########################################################################        

class Polymer:
    
    def __init__(self, idx, nbpp, r0, lbox, bid):
        """ INPUTS:
            idx: polymer index,
            nbpp: number of beads per polymer,
            r0: bond length
            lbox: box length
            bid: starting bead index in the polymers array"""
        
        ### store the beginning of the global bead index in the polymers array
        
        self.gid = bid
        
        ### atom and molecule types
        
        self.mol = np.ones((nbpp), dtype=np.int32)*idx
        self.tpe = np.ones((nbpp), dtype=np.int32)
        
        ### generate coordinates
        
        l = (nbpp-1)*r0                             # polymer length
        xi = l + random.random()*(lbox-2*l)         # random pos in box
        yi = l + random.random()*(lbox-2*l)
        ori = 2*np.pi*random.random()               # random orientation
        cosv = math.cos(ori)
        sinv = math.sin(ori)
        self.r = np.zeros((2, nbpp), dtype=np.float32)
        
        for j in xrange(nbpp):
            self.r[0,j] = xi + j*cosv
            self.r[1,j] = yi + j*sinv  

        ### generate bonds

        self.bonds = np.zeros((nbpp-1, 4), dtype=np.int32)
        for j in xrange(nbpp-1):
            self.bonds[j,0] = bid+j+1
            self.bonds[j,1] = 1
            self.bonds[j,2] = bid+j+1
            self.bonds[j,3] = bid+j+1+1 
        
        ### generate angles
        
        self.angles = np.zeros((nbpp-2, 5), dtype=np.int32)
        for j in xrange(nbpp-2):
            self.angles[j,0] = bid+j+1
            self.angles[j,1] = 1
            self.angles[j,2] = bid+j+1
            self.angles[j,3] = bid+j+1+1
            self.angles[j,4] = bid+j+1+2
        
        return

###########################################################################
        
def write_lammps_input(polymers, sim):
    """ write a lammps input file"""
    
    # open file for reading
    ofile = open('input.data', 'w') 
    # comment
    ofile.write('Ensemble of ' + str(nrod) + ' rods with ' + str(nbpr) + ' beads\n\n') 
    # number of atoms, bons, and angles
    natoms = len(x)
    nbonds = len(bonds)
    nangles = len(angles)
    ofile.write(str(natoms) + ' atoms\n')
    ofile.write(str(nbonds) + ' bonds\n')
    ofile.write(str(nangles) + ' angles\n\n')
    # atom, bond, and angle types
    #ofile.write(str(nbpr)+ ' atom types\n')
    ofile.write('1 atom types\n')
    ofile.write('1 bond types\n')
    ofile.write('1 angle types\n\n')
    # box dimensions
    ofile.write('0.0 ' + str(lbox) + ' xlo xhi\n')
    ofile.write('0.0 ' + str(lbox) + ' ylo yhi\n')
    ofile.write('-5.0 5.0 zlo zhi\n\n')
    # masses section
    ofile.write('Masses\n\n')
    ofile.write('1 1\n')
#    for i in range(nbpr):
#	ofile.write(str(i+1) + ' 1.0\n')
    ofile.write('\n')
    # atoms section
    ofile.write('Atoms\n\n')
    for i in range(natoms):
        # atom-ID, molecule-ID, atom-type, x,y,z
        ofile.write(str(i + 1) + ' ' + str(mol[i]) + ' ' + str(tpe[i]) +  ' ' + str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + '\n')
    ofile.write('\n')
    # bonds section
    ofile.write('Bonds\n\n')
    for i in range(nbonds):
        ofile.write(str(bonds[i,0]) + ' ' + str(bonds[i,1]) + ' ' + str(bonds[i,2]) + ' ' + str(bonds[i,3]) +  '\n')
    ofile.write('\n')
    # angles section
    ofile.write('Angles\n\n')
    for i in range(nangles):
        ofile.write(str(angles[i,0]) + ' ' + str(angles[i,1]) + ' ' + str(angles[i,2]) + ' ' + str(angles[i,3]) + ' ' + str(angles[i,4]) + '\n')
    ofile.write('\n')
    # close the file
    ofile.close()
    # print stats
    print 'Created box with',str(nrod),'rods'
    
    return
    
###########################################################################

def get_input_opts():
    """ get input options"""
 
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--density", \
                        type=float, help="Packing fraction")  
    parser.add_argument("-ll", "--llong", \
                        type=float, help="Length of the long polymers")
    parser.add_argument("-ls", "--lshort", 
                        type=float, help="Length of the short polymers")    
    parser.add_argument("-npl", "--npolsl", 
                        type=int, help="Number of long polymers in the simulation")        
    args = parser.parse_args()
    
    return args.llong, args.lshort, args.density, args.npolsl
    
###########################################################################

def main():
    
    ### specify parameters
    
    lpol_long, lpol_short, density, npols_long = get_input_opts()   
    sim = Simulation(lpol_long, lpol_short, density, npols_long)
    polymers = sim.gen_polymers()

    ### write the input file
    
    write_lammps_input(polymers, sim)
#    write_dummy_xyz(x, y, z)
    
    return

###########################################################################

if __name__ == '__main__':
    main()
    
###########################################################################    
