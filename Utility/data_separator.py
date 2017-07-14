#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 18:00:42 2017

@author: duman
"""

""" Separate data based on command line arguments,
    for cells, highdensefils and bidispersefils simulations"""

##############################################################################

import read_write
import misc_tools

##############################################################################

def gen_folder_path(base, sep, phaseparams):
    """ generate contiguous folder path address from the parameters, 
    sep should be the separator like _ or / """
    
    k = 0
    folder_path = ''
    for (key, value) in phaseparams:
        if k == 0:
            folder_path += base + key + '_' + str(value)
        else:
            folder_path += sep + key + '_' + str(value)
        k += 1
                
    return folder_path

##############################################################################
    
def gen_cell_data_address(eps, fp, areak, kappa, \
                          analysis, dbase, analysisdbase):
    """ data and analysis folders adress generator, 
    specific for cell work"""
    
    path1 = 'eps_' + str(eps) + \
        '/fp_' + str(fp) + '/areak_' + str(areak)  + \
        '/kappa_' + str(kappa)  
    path2 = analysis + '/' + analysis + \
        '_eps_' + str(eps) + '_fp_' + str(fp) + \
        '_areak_' + str(areak) + '_kappa_' + str(kappa) + '.txt'          
    datafolder = dbase + path1 + '/'
    analysisfile = analysisdbase + path2  

    return datafolder, analysisfile

##############################################################################
    
def gen_fil_data_address(density, kappa, fp, \
                          analysis, dbase, analysisdbase):
    """ data and analysis folders adress generator, 
    specific for filament work"""
    
    path1 = 'density_' + str(density) + '/kappa_' + str(kappa) + \
        '/fp_' + str(fp)
    path2 = analysis + '/' + analysis + '_density_' + str(density) + \
        '_kappa_' + str(kappa) + '_fp_' + str(fp) + '.txt'          
    datafolder = dbase + path1 + '/'
    analysisfile = analysisdbase + path2  

    return datafolder, analysisfile
    
##############################################################################

def read_cell_param_data(args):
    """ read cell parameter data"""

    legend_param = {}
    fixed_params = {}
    if args.kappa == -1:
        legend_param["kappa"] = [1.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0]
        fixed_params["eps"] = args.eps
        fixed_params["areak"] = args.areak
        fixed_params["fp"] = args.fp            
    elif args.fp == -1:
        legend_param["fp"] = [0.5, 1.0, 3.0, 5.0]  
        fixed_params["eps"] = args.eps
        fixed_params["areak"] = args.areak
        fixed_params["kappa"] = args.kappa   
    elif args.eps == -1:
        legend_param["eps"] = [0.05, 0.5, 1.0, 3.0, 5.0, 10.0, 20.0]
        fixed_params["kappa"] = args.kappa
        fixed_params["areak"] = args.areak
        fixed_params["fp"] = args.fp   
    elif args.areak == -1:
        legend_param["areak"] = [1.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0]
        fixed_params["eps"] = args.eps
        fixed_params["kappa"] = args.kapa
        fixed_params["fp"] = args.fp   
    
    return legend_param, fixed_params

##############################################################################
   
def read_highdens_fil_param_data(args):
    """ read filament parameter data,
    specific for highdens filaments data"""

    legend_param = {}
    fixed_params = {}
    if args.kappa == -1:
        #legend_param["kappa"] = [10.0, 20.0, 40.0, 75.0, 150.0]
        legend_param["kappa"] = [20.0, 40.0, 80.0, 100.0]
        fixed_params["fp"] = args.fp
        fixed_params["density"] = args.density
    elif args.fp == -1:
        #legend_param["fp"] = [0.0004, 0.004, 0.01, 0.04]  
        legend_param["fp"] = [0.0001, 0.001, 0.0024, 0.01] 
        fixed_params["density"] = args.density
        fixed_params["kappa"] = args.kappa   
    elif args.density == -1:
        legend_param["density"] = [0.8]
        fixed_params["kappa"] = args.kappa
        fixed_params["fp"] = args.fp   
    
    return legend_param, fixed_params

##############################################################################
   
def read_bidisperse_fil_param_data(args):
    """ read filament parameter data,
    specific for bidisperse filaments data"""

    legend_param = {}
    fixed_params = {}
    if args.kappa == -1:
        legend_param["kappa"] = [2.5, 25.0, 400.0]
        fixed_params["fp"] = args.fp
        fixed_params["density"] = args.density
    elif args.fp == -1:
        legend_param["fp"] = [0.024, 0.24, 2.4]  
        fixed_params["density"] = args.density
        fixed_params["kappa"] = args.kappa   
    elif args.density == -1:
        legend_param["density"] = [0.8]
        fixed_params["kappa"] = args.kappa
        fixed_params["fp"] = args.fp   
    
    return legend_param, fixed_params
    
##############################################################################    
    
def rearrange_cell_data(args, legend_param, fixed_params, read_analysis_func):
    """ rearrange the data to make it ready for plotting,
    namely encapsulate the data in a dictionary
    based on the chosen parameters,
    specific for cell work"""

    data = {}       # carries the data per parameter set
    sims = {}       # carries the simulation information per parameter set

    param_choice = legend_param.keys()[0]
    for p in legend_param[param_choice]:
        
        if param_choice == 'areak':
            datafolder, analysisfile = gen_cell_data_address(args.eps, \
                args.fp, p, args.kappa, args.analysisname, \
                args.folder, args.analysisbase)
        elif param_choice == 'eps':
            datafolder, analysisfile = gen_cell_data_address(p, \
                args.fp, args.areak, args.kappa, args.analysisname, \
                args.folder, args.analysisbase)
        elif param_choice == 'fp':            
            datafolder, analysisfile = gen_cell_data_address(args.eps, \
                p, args.areak, args.kappa, args.analysisname, \
                args.folder, args.analysisbase)
        elif param_choice == 'kappa':            
            datafolder, analysisfile = gen_cell_data_address(args.eps, \
                args.fp, args.areak, p, args.analysisname, \
                args.folder, args.analysisbase)
            
        sims[p] = read_write.read_sim_info(datafolder)
        out = read_analysis_func(analysisfile)
        if type(out) == tuple:
            if len(out) == 2:
                x, y = out
        else:
            y = out
        data[p] = y

    if type(out) == tuple:
        return x, data, param_choice, sims   
    else:
        return legend_param[param_choice], data, param_choice, sims    
    
##############################################################################    
    
def rearrange_fil_data(args, legend_param, fixed_params, read_analysis_func):
    """ rearrange the data to make it ready for plotting,
    namely encapsulate the data in a dictionary
    based on the chosen parameters,
    specific for filament work"""

    data = {}       # carries the data per parameter set
    sims = {}       # carries the simulation information per parameter set
    xp = {}         # carries the x axis of the data per parameter set

    param_choice = legend_param.keys()[0]
    for p in legend_param[param_choice]:
        
        if param_choice == 'density':
            datafolder, analysisfile = gen_fil_data_address(p, \
                args.kappa, args.fp, args.analysisname, \
                args.folder, args.analysisbase)
        elif param_choice == 'kappa':
            datafolder, analysisfile = gen_fil_data_address(args.density, \
                p, args.fp, args.analysisname, \
                args.folder, args.analysisbase)
        elif param_choice == 'fp':            
            datafolder, analysisfile = gen_fil_data_address(args.density, \
                args.kappa, p, args.analysisname, \
                args.folder, args.analysisbase)
            
        sims[p] = read_write.read_sim_info(datafolder)
        out = read_analysis_func(analysisfile)
        if type(out) == tuple:
            if len(out) == 2:
                x, y = out
                xp[p] = x
        else:
            y = out
        data[p] = y

    if type(out) == tuple:
        return xp, data, param_choice, sims   
    else:
        return legend_param[param_choice], data, param_choice, sims    
    
##############################################################################    

def gen_labels_for_cells(param_choice, sim):
    """ generate labels according to the parameter choice for helping in plots,
    specific for cell work"""
    
    eps_sym = r'$\epsilon$'
    f_sym = r'$Pe$'
    ka_sym = r'$\kappa_A$'
    kb_sym = r'$\kappa_B$'
    
    if param_choice == 'eps':
        xlab = eps_sym
        label = f_sym + "=" + str(sim.pe) + "," + \
            ka_sym + "=" + str(sim.areak) + "," + \
            kb_sym + "=" + str(sim.kappa)
    elif param_choice == 'fp':
        xlab = f_sym
        label = eps_sym + "=" + str(sim.eps) + "," + \
            ka_sym + "=" + str(sim.areak) + "," + \
            kb_sym + "=" + str(sim.kappa)        
    elif param_choice == 'areak':
        xlab = ka_sym
        label = eps_sym + "=" + str(sim.eps) + "," + \
            f_sym + "=" + str(sim.pe) + "," + \
            kb_sym + "=" + str(sim.kappa)        
    elif param_choice == 'kappa':
        xlab = kb_sym
        label = eps_sym + "=" + str(sim.eps) + "," + \
            f_sym + "=" + str(sim.pe) + "," + \
            ka_sym + "=" + str(sim.areak)        
    else:
        raise ValueError("Parameter choice is non-existing!")
    
    return xlab, label

##############################################################################

def gen_labels_for_fils(param_choice, sim):
    """ generate labels according to the parameter choice for helping in plots,
    specific for filaments works"""
    
    dens_sym = r'$\phi$'
    f_sym = r'$Pe$'
    k_sym = r'$\xi/L$'
    
    if param_choice == 'density':
        xlab = dens_sym
        label = k_sym + "=" + str(sim.xil) + "," + \
           f_sym + "=" + str(sim.pe) 
    elif param_choice == 'fp':
        xlab = f_sym
        label = dens_sym + "=" + str(sim.density) + "," + \
            k_sym + "=" + str(sim.xil)              
    elif param_choice == 'kappa':
        xlab = k_sym
        label = dens_sym + "=" + str(sim.density) + "," + \
            f_sym + "=" + str(sim.pe)
    else:
        raise ValueError("Parameter choice is non-existing!")
    
    return xlab, label
    
##############################################################################
    