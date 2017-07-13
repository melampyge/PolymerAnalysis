#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 14:40:40 2017

@author: duman
"""

""" generate povray rendered images of timeframes,
    color code can be polymer id or bond orientation"""

##############################################################################

### example command line arguments: 
###    -fl=/local/duman/SIMULATIONS/Bidisperse_Filaments/.../ 
###         -sb=/usr/users/iff_th2/duman/Bidisperse_Filaments/IMAGES/
###             -ti=10 -tf=2000 -c=id

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import read_write
import misc_tools
import data_structures
import vapory
import matplotlib.cm as cm
import matplotlib.colors as mplcolors 

##############################################################################

def get_args():
    """ get the command line arguments"""
    
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
    args = parser.parse_args()
    
    return args
    
##############################################################################
    
def gen_img_settings_quality(l):
    """ generate the general povray settings"""
    
    lhalf = 0.5*l
    
    ### sphere radius
    
    sphere_radius = 0.7
    
    ### RESOLUTION
    
    img_widthpx = 1024
    img_heightpx = 1024

    ### includes and defaults

    povray_includes = ["colors.inc", "textures.inc", "shapes.inc"]
    povray_defaults = [vapory.Finish( 'ambient', 0.1, \
	     			  'diffuse', 0.65, \
		    		  'specular', 0.5, \
			    	  'shininess', 0.53, \
				  'opacity', 1.0)]


    ### light sources

    sun1 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', 'White')
    sun2 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', [0.7, 0.7, 0.7])

    ### background

    background = vapory.Background('color', [1,1,1])

    ### camera

    povray_cam = vapory.Camera('location', [lhalf, lhalf, -1.01*lhalf], \
                               'look_at', [lhalf,lhalf,0], 'angle', 90)

    ### text
    # If desired include this in the povray_objects - array declared in the loop
    #text1 = \
    #vapory.Text( 'ttf', '"timrom.ttf"' ,'"Division:"', 0.01, 0.0, \
    #'scale', [0.5,0.5,0.5],'rotate', \
    #[0,90,0], 'translate' , [0.0 , 15.0+2.75-1 , 15.0+1.5], \
    #vapory.Pigment('Black') ) 

    ### render quality

    quality = 10
    
    return sphere_radius, img_widthpx, img_heightpx, povray_includes, \
        povray_defaults, sun1, sun2, background, povray_cam, quality

##############################################################################

def gen_colors_based_on_id(nbeads, npols, pid):
    """ generate an array with colors for the spheres,
    based on the polymer identities"""

    # use colormap (blue, red, green)
    #  blue 0.0, 0.0, 1.0
    #  red: 1.0, 0.0, 0.0
    #  green: 0.0, 1.0, 1.0

    my_cmap = cm.get_cmap('jet')
    norm = mplcolors.Normalize(0, npols)
    
    sphere_rgbcolor = []
    for i in xrange(nbeads):
        idx = pid[i]
        si = my_cmap(norm(idx))
        sphere_rgbcolor.append([si[0], si[1], si[2]])

    return sphere_rgbcolor
    
##############################################################################

def gen_colors_based_on_orient(nbeads, npols, ori):
    """ generate an array with colors for the spheres,
    based on the bond orientations"""

    # use colormap (blue, red, green)
    #  blue 0.0, 0.0, 1.0
    #  red: 1.0, 0.0, 0.0
    #  green: 0.0, 1.0, 1.0

    my_cmap = cm.get_cmap('hsv')
    norm = mplcolors.Normalize(0, 2.*np.pi)
    
    sphere_rgbcolor = []
    for i in xrange(nbeads):
        idx = ori[i]
        si = my_cmap(norm(idx))
        sphere_rgbcolor.append([si[0], si[1], si[2]])

    return sphere_rgbcolor    
       
##############################################################################

def plot_frames(beads, sim, ti, tf, savebase, colorid):
    """ plot frames within the specified time window"""
    
    ### define the color for the spheres

    print 'defining colors'
    if colorid == "id":
        sphere_rgbcolor = gen_colors_based_on_id(sim.nbeads, sim.npols, beads.pid)
    elif colorid == "orient":
        sphere_rgbcolor = gen_colors_based_on_orient(sim.nbeads, sim.npols, beads.ori)        

    ### create povray settings

    print 'creating povray settings'
    sphere_radius, img_widthpx, img_heightpx, povray_includes, \
        povray_defaults, sun1, sun2, background, povray_cam, quality \
            = gen_img_settings_quality(sim.lx)
    
    zi = np.zeros((sim.nbeads), dtype=np.float32)
        
    ### set general plot properties

    os.system("mkdir -p " + savebase)
    savebase = misc_tools.gen_contiguous_folder_path(savebase, sim.density, \
                                                     sim.kappa, sim.fp)
    os.system("mkdir -p " + savebase)
    
    ### plot the frames
    
    for step in range(ti, tf):
        
        time = step*sim.dt
        print 'Step / Total : ', step, tf
        
        ### create povray items
        
        print 'generating povray item'
        particles = vapory.Object( \
            vapory.Union( \
                *[ vapory.Sphere([beads.xi[step, 0, j], beads.xi[step, 1, j],zi[j]], \
                    sphere_radius, vapory.Texture( \
                        vapory.Pigment('color', sphere_rgbcolor[j]), \
                            vapory.Finish('phong',1)) ) for j in range(0, sim.nbeads ) ] ) )

        ### generate povray objects

        print 'generating povray objects'
        povray_objects = [sun1, sun2, background, particles]
        ### create the scene
        scene = vapory.Scene( camera = povray_cam,
                       objects = povray_objects, 
                       included = povray_includes, 
                       defaults = povray_defaults )
                       
        ### render image
                           
        print 'rendering scene'
        savename = "pov-frame-" + "{0:05d}".format(int(step)) + ".png"
        scene.render(outfile=savename, width=img_widthpx, height=img_heightpx, \
            antialiasing=0.001, quality=quality, remove_temp=True)
            
        ### move the image to the correct destination
            
        os.system('mv ' + savename + ' ' + savebase)
        
    return
        
##############################################################################

def main():

    ### get the command line arguments
    
    args = get_args()
    
    ### read the data and general information from the folder
    
    beads, pols, sim = read_write.read_h5_file(args.folder)
    print "folder = ", args.folder
    
    print "Calculating image positions of beads"
    beads.xi = misc_tools.get_img_pos(beads.xu, sim.lx)

    if args.colorid == "orient":
        print "Calculating bond orientations"
        beads.calc_bond_orientations(beads.xu, \
                        sim.nsteps, sim.nbeads, sim.npols, sim.nbpp)
    elif args.colorid == "id":
        print "Generating polymer identities of beads"
        beads.get_pol_ids(sim.nbeads, sim.nbpp)
        
    ### plot the data in the given time window
    
    plot_frames(beads, sim, args.init_time, args.fin_time, \
                args.savebase, args.colorid)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
