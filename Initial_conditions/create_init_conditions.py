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
import os

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
                                        
        self.npols['short'] = int(self.nbeads['short']/self.nbeads_per_pol['short'])
        self.nbeads['short'] = self.npols['short']*self.nbeads_per_pol['short']

        self.totnbeads = self.nbeads['long'] + self.nbeads['short']
                                        # total number of beads        
        self.totnpols = self.npols['long'] + self.npols['short']
                                        # total number of polymers
        
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
            r0: bond length,
            lbox: box length,
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
            self.bonds[j,0] = j+1
            self.bonds[j,1] = 1
            self.bonds[j,2] = bid+j+1
            self.bonds[j,3] = bid+j+1+1 
        
        ### generate angles
        
        self.angles = np.zeros((nbpp-2, 5), dtype=np.int32)
        for j in xrange(nbpp-2):
            self.angles[j,0] = j+1
            self.angles[j,1] = 1
            self.angles[j,2] = bid+j+1
            self.angles[j,3] = bid+j+1+1
            self.angles[j,4] = bid+j+1+2
        
        return

###########################################################################
        
def write_lammps_input(npols, x, y, z, bonds, angles, mol, tpe, lbox, folder):
    """ write a lammps input file"""
    
    ### open file for writing
    
    ofile = open(folder + 'input.data', 'w') 
    
    ### comment
    
    ofile.write('Ensemble of ' + str(npols) + \
                ' polymers with bidispersity\n\n') 
    
    ### number of atoms, bonds, and angles
    
    natoms = len(x)
    nbonds = len(bonds)
    nangles = len(angles)
    ofile.write(str(natoms) + ' atoms\n')
    ofile.write(str(nbonds) + ' bonds\n')
    ofile.write(str(nangles) + ' angles\n\n')
    
    ### atom, bond, and angle types

    ofile.write('1 atom types\n')
    ofile.write('1 bond types\n')
    ofile.write('1 angle types\n\n')
    
    ### box dimensions
    ofile.write('0.0 ' + str(lbox) + ' xlo xhi\n')
    ofile.write('0.0 ' + str(lbox) + ' ylo yhi\n')
    ofile.write('-5.0 5.0 zlo zhi\n\n')
    
    ### masses section
    
    ofile.write('Masses\n\n')
    ofile.write('1 1\n')
    ofile.write('\n')
    
    ### atoms section
    
    ofile.write('Atoms\n\n')
    for i in range(natoms):
        
        ### atom-ID, molecule-ID, atom-type, x, y, z
        
        ofile.write(str(i+1) + ' ' + str(mol[i]) + ' ' + \
                    str(tpe[i]) +  ' ' + str(x[i]) + ' ' + \
                    str(y[i]) + ' ' + str(z[i]) + '\n')
    ofile.write('\n')
    
    ### bonds section
    
    ofile.write('Bonds\n\n')
    for i in range(nbonds):
        ofile.write(str(bonds[i,0]) + ' ' + \
                    str(bonds[i,1]) + ' ' + str(bonds[i,2]) + \
                    ' ' + str(bonds[i,3]) +  '\n')
    ofile.write('\n')
    
    ### angles section
    
    ofile.write('Angles\n\n')
    for i in range(nangles):
        ofile.write(str(angles[i,0]) + ' ' + \
                    str(angles[i,1]) + ' ' + str(angles[i,2]) + \
                    ' ' + str(angles[i,3]) + ' ' + \
                    str(angles[i,4]) + '\n')
    ofile.write('\n')
    
    ### close the file
    
    ofile.close()
    
    ### print stats
    
    print 'Created box with ', str(npols), ' polymers'
    
    return

###########################################################################

def write_dummy_xyz(x, y, z, folder):
    """ write a xyz file for testing purposes"""
    
    natoms = len(x)
    ofile = open(folder + 'dummy.xyz', 'w')
    ofile.write(str(natoms) + '\n')
    ofile.write('comment line\n')
    for i in range(natoms):
        ofile.write('C\t' + str(x[i]) + '\t' + str(y[i]) + '\t' + str(z[i]) + '\n')
    ofile.close()
    
    return
    
###########################################################################

def serialize_data(polymers, sim):
    """ reformat the data into arrays and primitive data types
    this is a quick fix, that needs to be changed!"""

    ### allocation
    
    x = np.zeros((sim.totnbeads), dtype=np.float32)
    y = np.zeros((sim.totnbeads), dtype=np.float32)
    z = np.zeros((sim.totnbeads), dtype=np.float32)
    mol = np.zeros((sim.totnbeads), dtype=np.int32)
    tpe = np.ones((sim.totnbeads), dtype=np.int32)
    
    ### coordinates
    
    k = 0
    for n in xrange(sim.totnpols):
        mol[k] = n
        if n < sim.npols['long']:
            nbeads_per_pol = sim.nbeads_per_pol['long']
        else:
            nbeads_per_pol = sim.nbeads_per_pol['short']            
        for j in xrange(nbeads_per_pol):
            x[k] = polymers[n].r[0, j]
            y[k] = polymers[n].r[1, j]
            k += 1
        
    ### bonds
    
    nbonds = sim.npols['long']*(sim.nbeads_per_pol['long']-1) + \
            sim.npols['short']*(sim.nbeads_per_pol['short']-1)
    bonds = np.zeros((nbonds, 4), dtype=np.int32)
    
    k = 0
    for n in xrange(sim.totnpols):
        for j in xrange(len(polymers[n].bonds[:, 0])):
            bonds[k+j, 0] = polymers[n].bonds[j, 0] + k
            bonds[k+j, 1] = polymers[n].bonds[j, 1]
            bonds[k+j, 2] = polymers[n].bonds[j, 2]
            bonds[k+j, 3] = polymers[n].bonds[j, 3] 
        k += len(polymers[n].bonds[:, 0])
    
    ### angles
    
    nangles = sim.npols['long']*(sim.nbeads_per_pol['long']-2) + \
            sim.npols['short']*(sim.nbeads_per_pol['short']-2)
    angles = np.zeros((nangles, 5), dtype=np.int32)
    
    k = 0
    for n in xrange(sim.totnpols):
        for j in xrange(len(polymers[n].angles[:, 0])):
            angles[k+j, 0] = polymers[n].angles[j, 0] + k
            angles[k+j, 1] = polymers[n].angles[j, 1]
            angles[k+j, 2] = polymers[n].angles[j, 2] 
            angles[k+j, 3] = polymers[n].angles[j, 3] 
            angles[k+j, 4] = polymers[n].angles[j, 4]
        k += len(polymers[n].angles[:, 0])
        
    return sim.totnpols, x, y, z, bonds, angles, mol, tpe, sim.box_length
    
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
    parser.add_argument("-fd", "--folder",
                        type=str, help="Specify the folder to save the data inside")      
    args = parser.parse_args()
    
    os.system("mkdir -p " + args.folder)
    
    return args.llong, args.lshort, args.density, args.npolsl, args.folder
    
###########################################################################

def main():
    
    ### specify parameters
    
    lpol_long, lpol_short, density, npols_long, folder = get_input_opts()   
    sim = Simulation(lpol_long, lpol_short, density, npols_long)
    polymers = sim.gen_polymers()

    ### write the input file
    
    npols, x, y, z, bonds, angles, mol, tpe, lbox = serialize_data(polymers, sim)
    write_lammps_input(npols, x, y, z, bonds, angles, mol, tpe, lbox, folder)
    write_dummy_xyz(x, y, z, folder)
    
    return

###########################################################################

if __name__ == '__main__':
    main()
    
###########################################################################    
