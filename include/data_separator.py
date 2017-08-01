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
import data_structures
import os

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

class Separator:
    """ data separator for plotting purposes"""

    def __init__(self, args, read_fnc):

        self.simtype = args.simtype
        if self.simtype == "highdensfilaments":
            self.legend_param, self.fixed_param = self.read_highdens_fil_param_data(args)
        elif self.simtype == "bidispersefilaments":
            self.legend_param, self.fixed_param = self.read_bidisperse_fil_param_data(args)
        elif self.simtype == "cells":
            self.legend_param, self.fixed_param = self.read_cell_param_data(args)
        else:
            raise ValueError("The chosen simtype does not exist!")

        self.out = self.read_analysis_data(args, read_fnc)

        return

    ##############################################################################

    def read_cell_param_data(self, args):
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

    def read_highdens_fil_param_data(self, args):
        """ read filament parameter data,
        specific for highdens filaments data"""

        legend_param = {}
        fixed_params = {}
        if args.kappa == -1:
            legend_param["kappa"] = [6.0, 12.0, 24.0, 60.0]
            fixed_params["fp"] = args.fp
            fixed_params["density"] = args.density
        elif args.fp == -1:
            legend_param["fp"] = [0.001, 0.01, 0.027, 0.1]
            fixed_params["density"] = args.density
            fixed_params["kappa"] = args.kappa
        elif args.density == -1:
            legend_param["density"] = [0.8]
            fixed_params["kappa"] = args.kappa
            fixed_params["fp"] = args.fp

        return legend_param, fixed_params

    ##############################################################################

    def read_bidisperse_fil_param_data(self, args):
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

    def gen_cell_address(self, eps, fp, areak, kappa, \
                         analysis, dbase, analysisdbase):
        """ data and analysis folders address generator,
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

    def gen_fil_address(self, density, kappa, fp, \
                        analysis, dbase, analysisdbase):
        """ data and analysis folders address generator,
        specific for filament work"""

        path1 = 'density_' + str(density) + '/kappa_' + str(kappa) + \
            '/fp_' + str(fp)
        path2 = analysis + '/' + analysis + '_density_' + str(density) + \
            '_kappa_' + str(kappa) + '_fp_' + str(fp) + '.txt'
        datafolder = dbase + path1 + '/'
        analysisfile = analysisdbase + path2

        return datafolder, analysisfile

    ##############################################################################

    def read_analysis_data(self, args, read_analysis_func):
        """ read all the analysis data and
        rearrange the data to make it ready for plotting,
        namely encapsulate the data in a dictionary
        based on the chosen parameters"""

        data = {}       # carries the data per parameter set
        sims = {}       # carries the simulation information per parameter set
        xp = {}         # carries the x axis of the data per parameter set

        param_choice = self.legend_param.keys()[0]

        for p in self.legend_param[param_choice]:

            if self.simtype != "cells":
                if param_choice == 'density':
                    datafolder, analysisfile = self.gen_fil_address(p, \
                        args.kappa, args.fp, args.analysisname, \
                        args.folder, args.analysisbase)
                elif param_choice == 'kappa':
                    datafolder, analysisfile = self.gen_fil_address(args.density, \
                        p, args.fp, args.analysisname, \
                        args.folder, args.analysisbase)
                elif param_choice == 'fp':
                    datafolder, analysisfile = self.gen_fil_address(args.density, \
                        args.kappa, p, args.analysisname, \
                        args.folder, args.analysisbase)
                sims[p] = data_structures.SimulationFilaments(datafolder)
            else:
                if param_choice == 'eps':
                    datafolder, analysisfile = self.gen_cell_address(p, \
                        args.fp, args.areak, args.kappa, args.analysisname, \
                        args.folder, args.analysisbase)
                elif param_choice == 'kappa':
                    datafolder, analysisfile = self.gen_cell_address(args.eps, \
                        args.fp, args.areak, p, args.analysisname, \
                        args.folder, args.analysisbase)
                elif param_choice == 'fp':
                    datafolder, analysisfile = self.gen_cell_address(args.eps, \
                        p, ars.areak, args.kappa,  args.analysisname, \
                        args.folder, args.analysisbase)
                elif param_choice == 'areak':
                    datafolder, analysisfile = self.gen_cell_address(args.eps, \
                        args.fp, p, args.kappa, args.analysisname, \
                        args.folder, args.analysisbase)
                sims[p] = data_structures.SimulationCells(datafolder)

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

    def gen_labels_for_cells(self, param_choice, sim):
        """ generate labels according to the parameter choice for helping in plots,
        specific for cell work"""

        eps_sym = r'$\epsilon$'
        f_sym = r'$Pe$'
        ka_sym = r'$\kappa_A$'
        kb_sym = r'$\kappa_B$'
        eps_name = "epsilon"
        f_name = "fp"
        ka_name = "areak"
        kb_name = "kappa"

        if param_choice == 'eps':
            xlab = eps_sym
            label = f_sym + "=" + "{0:.2g}".format(sim.pe) + "," + \
                ka_sym + "=" + str(sim.areak) + "," + \
                kb_sym + "=" + str(sim.kappa)
            path = f_name + "=" + str(sim.fp) + "_" + \
                ka_name + "=" + str(sim.areak) + "_" + \
                kb_name + "=" + str(sim.kappa)
        elif param_choice == 'fp':
            xlab = f_sym
            label = eps_sym + "=" + str(sim.eps) + "," + \
                ka_sym + "=" + str(sim.areak) + "," + \
                kb_sym + "=" + str(sim.kappa)
            path = eps_name + "=" + str(sim.eps) + "_" + \
                ka_name + "=" + str(sim.areak) + "_" + \
                kb_name + "=" + str(sim.kappa)
        elif param_choice == 'areak':
            xlab = ka_sym
            label = eps_sym + "=" + str(sim.eps) + "," + \
                f_sym + "=" + "{0:.2g}".format(sim.pe) + "," + \
                kb_sym + "=" + str(sim.kappa)
            path = eps_name + "=" + str(sim.eps) + "_" + \
                f_name + "=" + str(sim.fp) + "_" + \
                kb_name + "=" + str(sim.kappa)
        elif param_choice == 'kappa':
            xlab = kb_sym
            label = eps_sym + "=" + str(sim.eps) + "," + \
                f_sym + "=" +"{0:.2g}".format(sim.pe)  + "," + \
                ka_sym + "=" + str(sim.areak)
            label = eps_name + "=" + str(sim.eps) + "_" + \
                f_name + "=" + str(sim.fp) + "_" + \
                ka_name + "=" + str(sim.areak)
        else:
            raise ValueError("Parameter choice is non-existing!")

        return xlab, label, path

    ##############################################################################

    def gen_labels_for_fils(self, param_choice, sim):
        """ generate labels according to the parameter choice for helping in plots,
        specific for filaments works"""

        dens_sym = r'$\phi$'
        f_sym = r'$Pe$'
        k_sym = r'$\xi/L$'
        dens_name = "density"
        f_name = "fp"
        k_name = "kappa"

        if param_choice == 'density':
            xlab = dens_sym
            label = k_sym + "=" + str(sim.xil) + "," + \
               f_sym + "=" + "{0:.2g}".format(sim.pe)
            path = k_name + "=" + str(sim.kappa) + "_" + \
                f_name + "=" + str(sim.fp)
        elif param_choice == 'fp':
            xlab = f_sym
            label = dens_sym + "=" + str(sim.density) + "," + \
                k_sym + "=" + str(sim.xil)
            path = dens_name + "=" + str(sim.density) + "_" + \
                k_name + "=" + str(sim.kappa)
        elif param_choice == 'kappa':
            xlab = k_sym
            label = dens_sym + "=" + str(sim.density) + "," + \
                f_sym + "=" + "{0:.2g}".format(sim.pe)
            path = dens_name + "=" + str(sim.density) + "_" + \
                f_name + "=" + str(sim.fp)
        else:
            raise ValueError("Parameter choice is non-existing!")

        return xlab, label, path

    ##############################################################################

