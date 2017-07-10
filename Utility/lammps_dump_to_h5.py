
""" Combine multiple dump files from different restart instances 
    into a single hdf5 file"""

### example command line arguments: 
###    -fd=/homea/ias2/duman/Bidisperse_Filaments/ -t=100000000
###         -d=0.8 -k=0.0 -f=0.0 
###         -dt=0.001 -ns=50000 -b=0.5 -s=1.0 
###         -ll=20.0 -ls=8.0 -npl=1000

##############################################################################

import argparse
import numpy as np
import os
import h5py
import subprocess
import data_structures
import misc_tools

##############################################################################

def read_contextual_info():
    """ read the contextual information provided by the user'
    specific for bidisperse polymers"""
    
    ### get the data folder and the last timestep info
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fd", "--folder", type=str, \
                        help="Folder containing data, as in /homea/ias2/duman/Bidisperse_Filaments/")
    parser.add_argument("-t", "--last_tstep", nargs="?", const="100000000", \
                            type=int, help="The last time step that is being searched for")
    parser.add_argument("-d", "--density", type=float, \
                        help="Packing fraction of the system")  
    parser.add_argument("-k", "--kappa", type=float, \
                        help="Bending rigidity")
    parser.add_argument("-f", "--fp", type=float, \
                        help="Propulsion force")
    parser.add_argument("-dt", "--timestep", type=float, \
                        help="Timestep of the simulation")
    parser.add_argument("-ns", "--nsamp", type=int, \
                        help="Sampling rate of data")
    parser.add_argument("-b", "--bl", type=float, \
                        help="Bond length of the simulation")    
    parser.add_argument("-s", "--sigma", type=float, \
                        help="Lennard Jones length")     
    parser.add_argument("-ll", "--llong", \
                        type=float, help="Length of the long polymers")
    parser.add_argument("-ls", "--lshort", \
                        type=float, help="Length of the short polymers")    
    parser.add_argument("-npl", "--npolsl", \
                        type=int, help="Number of long polymers in the simulation") 
    args = parser.parse_args()
    
    ### generate folder path
    
    folder = misc_tools.gen_folder_path(args.folder, args.density, args.kappa, args.fp)
    fl = folder + 'out1.dump'
    assert os.path.exists(fl), "\nOUT1.DUMP DOES NOT EXIST FOR: " + folder 

    print "Creating the compressed data format for the following file : " + fl

    ### determine number of beads per polymer
    
    args.nbpp_long = int(args.llong/args.bl+1)
    args.nbpp_short = int(args.lshort/args.bl+1) 
    
    ### determine number of beads per polymer type in total
    
    args.nbeads_long = args.npolsl*args.nbpp_long
    nbeads_short = args.nbeads_long
    
    ### determine number of short polymers
    
    args.npolss = int(nbeads_short/args.nbpp_short)
    
    ### determine number of beads per short polymers
    
    args.nbeads_short = args.npolss*args.nbpp_short
    
    ### determine total number of beads and polymers
    
    args.totnbeads = args.nbeads_long + args.nbeads_short
    args.totnpols = args.npolsl + args.npolss
        
    ### determine box information 
    
    area_of_pol = args.llong*args.sigma*args.npolsl + \
        args.lshort*args.sigma*args.npolss
    area_box = area_of_pol/args.density
    args.lx = np.sqrt(area_box)
    args.ly = np.sqrt(area_box)
    
    ### total number of steps
    
    args.nsteps = args.last_tstep/args.nsamp
    args.nsteps += 1
    
    ### restructure the data 
    
    nbpp = np.ones((args.totnpols), dtype=np.int32)
    nbpp[0:args.npolsl] *= args.nbpp_long
    nbpp[args.npolsl:] *= args.nbpp_short
    
    sim = data_structures.SimulationBidispersePolymers(folder, \
                                     args.timestep*args.nsamp, args.density, \
                                     args.nsteps, args.totnbeads, \
                                     args.totnpols, nbpp, args.bl, args.sigma, \
                                     args.lx, args.ly, args.kappa, args.fp)
    
    return sim, args.last_tstep
        
##############################################################################
    
def get_number_of_snaps(f, nbeads):
    """ determine the number of snapshots in the file --old version"""
    
    os.system('wc -l ' + f + ' > tmp.txt')
    ifile = open('tmp.txt')
    line = ifile.readline()
    line = line.split()
    nlines = int(line[0])
    nsnaps = nlines/(nbeads+9)
    ifile.close()
    os.system('rm tmp.txt')
    
    return nsnaps

##############################################################################
    
def get_number_of_snaps_alternative_trial(f, nbeads):
    """ determine the number of snapshots in the file"""
    
    ### THIS DOES NOT WORK YET!!
    
    nlines = int(subprocess.check_output(['wc', '-l', f]))
    print 'subprocess worked'
    print 'nlines = ', str(nlines)
    nsnaps = nlines/(nbeads+9)
    
    return nsnaps

##############################################################################
    
def read_pos(fl, x, mid, T, nbeads, nsnaps, lx, ly, checked, tstep_cnt):
    """ read the position data from the file and return the last timestep at the end"""

    ### read the positions unique per each tstep in a single dump file

    already_checked = False
    for snap in range(nsnaps):
        
        ### read the headers to check the uniqueness of current tstep
        
        fl.readline()
        line = fl.readline()
        line = line.split()
        tstep = int(line[0])
        
        ### finish if the last tstep is exceeded already
        
        if tstep > T:
            return tstep_cnt, tstep
        
        ### make sure the current tstep is unique
        
        if tstep not in checked:
            checked.append(tstep)
            tstep_cnt += 1
            already_checked = False
        else:
            print "ALREADY CHECKED " + str(tstep)             
            already_checked = True
    
        ### read the remaining part of the header part
        
        for j in range(7):
            fl.readline()
         
        ### read the positions per bead if the tstep is unique
        
        for j in range(nbeads):
            line = fl.readline()
            if already_checked:
                continue
            line = line.split()
            if len(line) < 2:
                print line
                continue
            bid = int(line[0]) - 1
            mid[bid] = int(line[1]) - 1
            x[tstep_cnt, 0, bid] = float(line[2])
            x[tstep_cnt, 1, bid] = float(line[3])

        print tstep_cnt, tstep

    return tstep_cnt, tstep

##############################################################################
    
def read_pos_from_dump_files(sim, T):
    """ read the position data of each dump file until the last tstep is reached"""

    ### generate file path and the total number of snapshots in the file
    
    current_file_number = 0
    x = np.zeros((sim.nsteps, 2, sim.nbeads), dtype=np.float32)
    mid = np.zeros((sim.nbeads), dtype=np.int32)
    tstep_cnt = -1
    checked = []
    tstep = 0
    
    while tstep < T:
        
        current_file_number += 1
        fpath = sim.folder + '/out' + str(current_file_number) + '.dump'
        assert os.path.exists(fpath), "out dump file does NOT exist for: " + fpath
        fl = open(fpath, 'r')
        nsnaps = get_number_of_snaps(fpath, sim.nbeads)
        
        ### read the positions unique per each tstep in a single dump file
        
        tstep_cnt, tstep = read_pos(fl, x, mid, T, sim.nbeads, nsnaps, sim.lx, sim.ly, \
                                    checked, tstep_cnt)
        fl.close()
        
    return x, mid
    
##############################################################################
    
def write_h5_file(xu, comu, sim):
    """ write data to hdf5 file,
    specific to bidisperse filament simulations"""
    
    ### file path
    
    fpath = sim.folder + 'out.h5'
    print "Writing data to ", fpath
    fl = h5py.File(fpath, 'w')
    
    ### positions of beads
    
    bead = fl.create_group('beads')
    bead.create_dataset('xu', (sim.nsteps, 2, sim.nbeads), data=xu, \
                        dtype=np.float32, compression='gzip') 
    
    ### polymer information
    
    pols = fl.create_group('pols')
    pols.create_dataset('comu', (sim.nsteps, 2, sim.npols), data=comu, \
                        dtype=np.float32, compression='gzip') 
    
    ### general information about the simulation
    
    info = fl.create_group('sim')
    info.create_dataset('nsteps', data=sim.nsteps)
    info.create_dataset('nbeads', data=sim.nbeads)
    info.create_dataset('npols', data=sim.npols)
    info.create_dataset('nbpp', data=sim.nbpp)
    
    ### simulation parameters  --- specific to analysis
    
    param = fl.create_group('params')
    param.create_dataset('density', data=sim.density)
    param.create_dataset('kappa', data=sim.kappa)
    param.create_dataset('fp', data=sim.fp)
    param.create_dataset('bl', data=sim.bl)
    param.create_dataset('sigma', data=sim.sigma)
    param.create_dataset('lx', data=sim.lx)
    param.create_dataset('ly', data=sim.ly)
    param.create_dataset('dt', data=sim.dt)
    
    fl.close()
    
    return

##############################################################################
    
def calculate_com_of_pols(xu, nsteps, nbpp, npols):   
    """ calculate the center of mass of polymers"""
    
    com = np.zeros((nsteps, 2, npols), dtype=np.float32)
    
    k = 0
    for j in range(npols):
        com[:, :, j] = np.mean(xu[:, :, k:k+nbpp[j]], axis=2)
        k += nbpp[j]
    
    return com

##############################################################################
    
def calculate_com_of_pols_one_liner(xu, nsteps, nbpp, nbeads):   
    """ calculate the center of mass of polymers"""
    
    splitted = np.split(xu, np.cumsum(nbpp)[:-1], axis=2)
    r = np.array([np.mean(sfil, axis=2) for sfil in splitted])
    com = np.swapaxes(np.swapaxes(r, 0, 1), 1, 2)
    
    return com
    
##############################################################################
    
def main():

    ### read information given by user about the simulation
    
    sim, last_tstep = read_contextual_info()
    
    ### read the bead data
    
    xu, mid = read_pos_from_dump_files(sim, last_tstep)
    
    ### generate polymer data
    
    com = calculate_com_of_pols_one_liner(xu, sim.nsteps, sim.nbpp, sim.nbeads)
    
    ### write the data in hdf5 format 
    
    write_h5_file(xu, com, sim)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################