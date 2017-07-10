#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 16:46:21 2017

@author: duman
"""

""" Plot snapshots of timeframes from the simulation"""

##############################################################################

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import h5py
import os
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import matplotlib.colors as mplcolors 

import read_write
import misc_tools
import data_structures
import plotter

import seaborn as sns
sns.set(style="white",context='paper',
        font_scale=1.2,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 30, 
            'font.family': 'sans',"figure.dpi":300,
            "xtick.major.size": 8, "ytick.major.size": 8,
            'grid.linestyle': '--'})   

##############################################################################

def set_color_pallette(beads, sim, colorid):
    """ set the color pallette"""
    
    quant_steps = 2056
    if colorid == "id":
        minval = 0.
        maxval = 2.*np.pi
        cmap_ax = plt.cm.get_cmap('jet', quant_steps)
        norm_ax = mpl.colors.Normalize(vmin=minval, vmax=maxval)      
        cidx = beads.pid
    elif colorid == "orient":
        minval = 0.
        maxval = 2.*np.pi
        cmap_ax = plt.cm.get_cmap('hsv', quant_steps)        
        norm_ax = mpl.colors.Normalize(vmin=0, vmax=2.*np.pi)  
        cidx = beads.ori
    
    return cidx, minval, maxval, cmap_ax, norm_ax
  
##############################################################################
    
def set_plot_props(sim):
    """ set plot properties"""
    
    full_box_downlim = -1.
    full_box_uplim = sim.lx+1.
    
    return full_box_downlim, full_box_uplim
    
##############################################################################

def plot_frames(beads, sim, ti, tf, savebase, colorid, savepdf):
    """ plot the selected timeframes"""
    
    ax0 = plotter.GeneralPlot()
    cid, minval, maxval, cmap_ax0, norm_ax0 = set_color_pallette(beads, \
                                              sim, colorid)
    full_box_downlim, full_box_uplim = set_plot_props(sim)
    os.system("mkdir -p " + savebase)
    savebase = misc_tools.gen_folder_path(savebase, sim.density, \
                                          sim.kappa, sim.fp)
    os.system("mkdir -p " + savebase)
    
    for step in np.arange(ti, tf):
        line0 = ax0.scatter(beads.xi[step, 0, :]/sim.bond_length, \
                            beads.xi[step, 1, :]/sim.bond_length, \
                            s=6.0, \
                            c=cidx, cmap=cmap_ax0, \
                            edgecolors='None', alpha=0.8, \
                            vmin=minval, vmax=maxval, \
                            norm=norm_ax0, rasterized=True) 
        ax0.set_axis('scaled')
        savepath = savebase + "frame-" + "{0:05d}".format(int(step))
        
        ### labels
    
        ax0.set_xlabel(r"$x/r_0$", fontsize=40)
        ax0.set_ylabel(r"$y/r_0$", fontsize=40)
    
        ### limits
    
        ax0.set_xlim((full_box_downlim, full_box_uplim))
        ax0.set_ylim((full_box_downlim, full_box_uplim))
        
        ### ticks
        
        #ax0.xaxis.set_ticks(full_box_ticks)
        #ax0.yaxis.set_ticks(full_box_ticks)
        ax0.tick_params(axis='both', which='major', labelsize=30)    
                
#        if colorid == 'orient':
#        
#            cax0 = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg+ax_len/3, \
#                             ax_len/4.6, ax_len/4.6], projection='polar')
#            xval = np.arange(-np.pi, np.pi, 0.01)
#            yval = np.ones_like(xval)
#            cax0.scatter(xval, yval, c=xval, s=300, \
#               cmap=plt.cm.get_cmap('hsv',quant_steps), norm=norm, linewidths=0)
#            cax0.set_xticks([])
#            cax0.set_yticks([])
#            cax0.set_title('$\\phi$',fontsize=20)
#            cax0.set_rlim([-1,1])
#            cax0.set_axis_off()
            
        ### save 
    
        if savepdf:
            plt.savefig(savepath+".pdf", dpi=300, bbox_inches='tight', pad_inches=0.08)
        else:            
            plt.savefig(savepath+".png", dpi=300, bbox_inches='tight', pad_inches=0.08)
            
        ax0.fig.clf()
        
    return

##############################################################################
    
def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fd", "--folder", \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Bidisperse_Filaments/.../")
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Bidisperse_Filaments/IMAGES/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Bidisperse_Filaments/IMAGES/")     
    parser.add_argument("-ti","--init_time", nargs="?", const=100, type=int, \
                        help="First frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=1900, type=int, \
                        help="Last frame of the video (in terms of frame number), you can also leave it empty")
    parser.add_argument("-c","--colorid", type=str, \
                        help ="Decide on the coloring -id or orient-")   
    parser.add_argument("-s","--savepdf", action="store_true", \
                        help ="Decide on the coloring -id or orient-")                                   
    args = parser.parse_args()
    
    ### read the data and general information from the folder
    
    beads, pols, sim = read_write.read_h5_file(args.folder, args.bending)
    print "folder = ", args.folder
    
    print "Calculating image positions of beads"
    beads.xi = misc_tools.get_img_pos(beads.xu, sim.lx)

    if args.colorid == "orient":
        print "Calculating bond orientations"
        beads.ori = misc_tools.calc_bond_orientations(beads.xu, \
                        sim.nsteps, sim.nbeads, sim.npols, sim.nbpp)
    elif args.colorid == "id":
        print "Generating polymer identities of beads"
        beads.pid = misc_tools.get_pol_ids(sim.nbeads, sim.nbpp)
        
    ### plot the data in the given time window
    
    plot_frames(beads, sim, args.init_time, args.fin_time, \
                args.savebase, args.colorid, args.savepdf)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################

