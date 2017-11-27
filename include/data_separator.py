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
from collections import OrderedDict
import numpy as np

##############################################################################

def assign_physicalvalues(params, sim):
    """ assign physical parameter values to fixed params
    based on a single simulation instance"""

    for param in params:
        param.assign_physicalvalue(sim)

    return

##############################################################################

def gen_path(params):
    """ generate path for plots
    titles are fixed_param
    NOTE THAT params should already be ordered"""

    k = 0
    path = ''
    for param in params:
        if k == 0:
            path += param.physicalparamname + '=' + param.physicalvalue
        else:
            path += "_" + param.physicalparamname + '=' + param.physicalvalue
        k += 1

    return path

##############################################################################

def gen_title(params):
    """ generate title for plots
    titles are fixed_param
    NOTE THAT params should already be ordered"""

    k = 0
    title = ''
    for param in params:
        if k == 0:
            title += param.physicalparamsymbol + '=' + param.physicalvalue
        else:
            title += ", " + param.physicalparamsymbol + '=' + param.physicalvalue
        k += 1

    return title

##############################################################################

def gen_label(param):
    """ generate a label-type path"""

    return param.physicalparamsymbol + "=" + param.physicalvalue

##############################################################################

def gen_tag(name, value, sep):
    """ generate a tag-type path"""

    return name + sep + str(value)

##############################################################################

class BendingRigidity:

    def __init__(self, simtype):

        self.simfoldername = 'kappa'
        self.name = 'kappa'
        self.physicalparamname = 'xil'
        self.symbol = r'$\kappa_B$'
        if simtype == "bifils" or simtype == "highfils":
            self.physicalparamsymbol = r'$\xi_p/L$'
        elif simtype == "cells":
            self.physicalparamsymbol = r'$\xi_p/R$'
        else:
            raise ValueError("Simulation type is chosen wrongly")

        return

    def assign_value(self, single_value):

        self.value = single_value

        return

    def assign_physicalvalue(self, sim):

        self.physicalvalue = "{0:.2f}".format(sim.xil)

        return

    def assign_values(self, list_of_values):

        self.values = list_of_values

        return

    def set_predef_values(self, simtype):

        if simtype == "cells":
            self.assign_values([1.0, 10.0, 50.0, 100.0, 500.0, 1000.0])
        elif simtype == "bifils":
            pass
        elif simtype == "highfils":
            self.assign_values([6.0, 12.0, 24.0, 60.0])
        else:
            raise ValueError("Simulation type is chosen wrongly")

        return

    def assign_order(self, order):

        self.order = order

        return

#############################################################################

class PropulsionForce:

    def __init__(self):

        self.simfoldername = 'fp'
        self.name = 'fp'
        self.physicalparamname = 'pe'
        self.symbol = r'$\f_{p}$'
        self.physicalparamsymbol = r'$Pe$'

        return

    def assign_value(self, single_value):

        self.value = single_value

        return

    def assign_physicalvalue(self, sim):

        self.physicalvalue = "{0:.2f}".format(sim.pe)

        return

    def assign_values(self, list_of_values):

        self.values = list_of_values

        return

    def set_predef_values(self, simtype):

        if simtype == "cells":
            self.assign_values([0.5, 1.0, 3.0, 5.0, 10.0, 20.0])
        elif simtype == "bifils":
            pass
        elif simtype == "highfils":
            self.assign_values([0.001, 0.01, 0.027, 0.1])
        else:
            raise ValueError("Simulation type is chosen wrongly")

        return

    def assign_order(self, order):

        self.order = order

        return

#############################################################################

class Density:

    def __init__(self):

        self.simfoldername = 'density'
        self.name = 'density'
        self.physicalparamname = 'density'
        self.symbol = r'$\phi$'
        self.physicalparamsymbol = r'$\phi$'

        return

    def assign_value(self, single_value):

        self.value = single_value

        return

    def assign_physicalvalue(self, sim):

        self.physicalvalue = "{0:.2f}".format(sim.density)

        return

    def assign_values(self, list_of_values):

        self.values = list_of_values

        return

    def set_predef_values(self, simtype):

        if simtype == "cells":
            self.assign_values([0.8])
        elif simtype == "bifils":
            self.assign_values([0.8])
        elif simtype == "highfils":
            self.assign_values([0.8])
        else:
            raise ValueError("Simulation type is chosen wrongly")

        return

    def assign_order(self, order):

        self.order = order

        return

#############################################################################

class AreaCompressionModulus:

    def __init__(self):

        self.simfoldername = 'areak'
        self.name = 'areak'
        self.physicalparamname = 'areak'
        self.symbol = r'$\kappa_A$'
        self.physicalparamsymbol = r'$\kappa_A$'

        return

    def assign_value(self, single_value):

        self.value = single_value

        return

    def assign_physicalvalue(self, sim):

        self.physicalvalue = "{0:.1f}".format(sim.areak)

        return

    def assign_values(self, list_of_values):

        self.values = list_of_values

        return

    def set_predef_values(self, simtype):

        if simtype == "cells":
            self.assign_values([1.0, 10.0, 50.0, 100.0, 500.0, 1000.0])
        else:
            raise ValueError("Simulation type is chosen wrongly")

        return

    def assign_order(self, order):

        self.order = order

        return

#############################################################################

class InterEnergy:

    def __init__(self):

        self.simfoldername = 'eps'
        self.name = 'eps'
        self.physicalparamname = 'eps'
        self.symbol = r'$\epsilon$'
        self.physicalparamsymbol = r'$\epsilon/k_{B}T$'

        return

    def assign_value(self, single_value):

        self.value = single_value

        return

    def assign_physicalvalue(self, sim):

        self.physicalvalue = "{0:.2f}".format(sim.eps/sim.kT)

        return

    def assign_values(self, list_of_values):

        self.values = list_of_values

        return

    def set_predef_values(self, simtype):

        if simtype == "cells":
            self.assign_values([0.05, 0.5, 1.0, 3.0, 5.0, 10.0, 20.0])
        else:
            raise ValueError("Simulation type is chosen wrongly")

        return

    def assign_order(self, order):

        self.order = order

        return

#############################################################################

def paramkeys_to_objects(key, simtype):

    if key == "eps":
        return InterEnergy()
    elif key == "fp":
        return PropulsionForce()
    elif key == "kappa":
        return BendingRigidity(simtype)
    elif key == "areak":
        return AreaCompressionModulus()
    elif key == "density":
        return Density()
    else:
        raise ValueError("Wrong key is given")
    return

#############################################################################

def paramkeys_to_values(obj, key):

    if key == "eps":
        return obj.eps
    elif key == "fp":
        return obj.fp
    elif key == "kappa":
        return obj.kappa
    elif key == "areak":
        return obj.areak
    elif key == "density":
        return obj.density
    else:
        raise ValueError("Wrong key is given")
    return

#############################################################################

def all_param_combinations(simtype):

    ### NOTE THAT this list should be ordered

    if simtype == "cells":
        return ["eps", "fp", "areak", "kappa"]
    elif simtype == "bifils":
        return ["density", "kappa", "fp"]
    elif simtype == "highfils":
        return ["density", "kappa", "fp"]
    else:
        raise ValueError("Simulation type is chosen wrongly")
    return

#############################################################################

class Separator:
    """ data separator and holder to facilitate
    access to analysis data for plotting purposes"""

    def __init__(self, args, read_fnc):

        self.sim_type = args.sim_type
        self.plot_type = args.plot_type
        self.set_param_context(args, read_fnc)

        return

    def set_param_context(self, args, read_fnc):
        """ set the chosen parameters
        in context with metadata"""

        ### set all the parameter name combinations
        # according to simulation type

        self.allparams = all_param_combinations(self.sim_type)

        ### set all the specialized parameter names
        # according to plot and simulation types

        self.fixed_param = []
        if self.plot_type == "2D":
            self.legend_param = paramkeys_to_objects(args.legend, self.sim_type)
            self.legend_param.set_predef_values(self.sim_type)
            self.control_param = None
            for param in self.allparams:
                if param != args.legend:
                    self.fixed_param.append(paramkeys_to_objects(param, self.sim_type))

        elif self.plot_type == "1D":
            self.legend_param = paramkeys_to_objects(args.legend, self.sim_type)
            self.legend_param.set_predef_values(self.sim_type)
            self.control_param = paramkeys_to_objects(args.control, self.sim_type)
            self.control_param.set_predef_values(self.sim_type)
            for param in self.allparams:
                if param != args.legend and param != args.control:
                    self.fixed_param.append(paramkeys_to_objects(param, self.sim_type))

        elif self.plot_type == "phase":
            self.legend_param = None
            self.control_param = []
            self.control_param.append(paramkeys_to_objects(args.control[0], self.sim_type))
            self.control_param[0].set_predef_values(self.sim_type)
            self.control_param.append(paramkeys_to_objects(args.control[1], self.sim_type))
            self.control_param[1].set_predef_values(self.sim_type)
            for param in self.allparams:
                if param not in args.control:
                    self.fixed_param.append(paramkeys_to_objects(param, self.sim_type))

        elif self.plot_type == "single":
            self.legend_param = None
            self.control_param = None
            for param in self.allparams:
                self.fixed_param.append(paramkeys_to_objects(param, self.sim_type))

        else:
            raise ValueError("Plot type is chosen wrongly")

        ### set ordered parameter list

        self.order_params()

        ### set actual values to parameters based on simulation type
        # and earlier choices

        self.set_param_values(args)

        ### set the bases for folders to access simulation
        # and analysis folders

        self.set_folder_context(args)

        ### read the analysis data and simulation data
        # according to the user choices

        self.read_analysis_data(read_fnc)

        ### set labels, titles, texts, etc
        # for plotting

        #self.set_plot_context()

        return

    def set_param_values(self, args):
        """ set the parameter values
        based on the generated context"""

        for param in self.fixed_param:
            value = paramkeys_to_values(args, param.name)
            param.assign_value(value)

        return

    def add_param_to_ordered_params(self, param, param_check, k):
        """ check and if it is at the correct place add parameter
        to the ordered list of parameters"""

        update = False
        if param == param_check.name:
            self.ordered_params.append(param_check)
            k += 1
            update = True

        return k, update

    def order_params(self):
        """ order the parameters according to the order
        as set in allparams"""

        k = -1
        self.ordered_params = []
        for param in self.allparams:
            if self.legend_param is not None:
                k, u = self.add_param_to_ordered_params(param, self.legend_param, k)
                if u:
                    self.legend_param.assign_order(k)
            if self.control_param is not None:
                if type(self.control_param) == list:
                    for j, controlparam in enumerate(self.control_param):
                        k, u = self.add_param_to_ordered_params(param, self.controlparam[j], k)
                        if u:
                            self.control_param[j].assign_order(k)
                else:
                    if param == self.control_param.name:
                        k, u = self.add_param_to_ordered_params(param, self.control_param, k)
                        if u:
                            self.control_param.assign_order(k)
            for j, fixedparam in enumerate(self.fixed_param):
                k, u = self.add_param_to_ordered_params(param, fixedparam, k)
                if u:
                    self.fixed_param[j].assign_order(k)

        return

    def set_folder_context(self, args):
        """ set the folder context
        for data access"""

        ### NOTE :
        # I should make these options set available
        # by command line options

        path_416 = '/local/duman/SIMULATIONS/'
        path_home = '/usr/users/iff_th2/duman/'
        self.simfolderbase = path_416
        self.datafolderbase = path_home
        self.savefolderbase = path_home

        if self.sim_type == "cells":
            self.simfolderbase += "Cells_in_LAMMPS/density_0.8/"
            self.datafolderbase += "Cells_in_LAMMPS/DATA/"
            self.savefolderbase += "Cells_in_LAMMPS/PLOTS/"

        elif self.sim_type == "bifils":
            self.simfolderbase += "Bidisperse_Filaments/Simulations/"
            self.datafolderbase += "Bidisperse_Filaments/DATA/"
            self.savefolderbase += "Bidisperse_Filaments/PLOTS/"

        elif self.sim_type == "highfils":
            self.simfolderbase += "HighDens_Filaments/Simulations/"
            self.datafolderbase += "HighDens_Filaments/DATA/"
            self.savefolderbase += "HighDens_Filaments/PLOTS/"

        else:
            raise ValueError("Simulation type is chosen wrongly")

        self.datafolderbase += args.analysis_name + "/" \
            + args.analysis_name + "_"
        self.savefolderbase += args.analysis_name + '/'

        return

    def read_analysis_data(self, read_fnc):
        """ read and parse analysis data"""

        ### set keys to hash the data
        # and set data as values to the keys

        if self.plot_type == "2D":
            self.key = self.legend_param
            self.read_data_for_one_key(read_fnc)

        elif self.plot_type == "1D":
            self.key = (self.legend_param, self.control_param)
            self.read_data_for_two_keys(read_fnc)

        elif self.plot_type == "phase":
            self.key = (self.control_param[0], self.control_param[1])
            self.read_data_for_two_keys(read_fnc)

        elif self.plot_type == "single":
            self.key = None
            self.read_data_single_instance(read_fnc)

        else:
            raise ValueError("Plot type is chosen wrongly")

        return

    def read_data_for_one_key(self, read_fnc):
        """ read data specific for single keyed hashing
        --2D plot type"""

        self.data = {}
        for j, val in enumerate(self.key.values):
            simfolder = self.simfolderbase
            datafolder = self.datafolderbase
            k = 0
            for param in self.ordered_params:
                if param.order == k:
                    if param != self.key:
                        simfolder += gen_tag(param.name, param.value, "_") \
                            + "/"
                        datafolder += gen_tag(param.name, param.value, "_")
                    else:
                        simfolder += gen_tag(param.name, param.values[j], "_") \
                            + "/"
                        datafolder += gen_tag(param.name, param.values[j], "_")

                    if k != len(self.ordered_params)-1:
                        datafolder += "_"
                    k += 1

            if self.sim_type == "cells":
                sim = data_structures.SimulationCells(simfolder)
            else:
                sim = data_structures.SimulationFilaments(simfolder)

            datafile = datafolder + ".txt"
            analysisdata = read_fnc(datafile)
            self.data[val] = Analysis(sim, analysisdata)

        return

    def read_data_for_two_keys(self, read_fnc):
        """ read data specific for two keyed hashing
        --1D and phase plot types"""

        self.data = {}
        for i, vali in enumerate(self.key[0].values):
            for j, valj in enumerate(self.key[1].values):
                simfolder = self.simfolderbase
                datafolder = self.datafolderbase
                k = 0
                for param in self.ordered_params:
                    if param.order == k:
                        if param not in self.key:
                            simfolder += gen_tag(param.name, param.value, "_") \
                                + "/"
                            datafolder += gen_tag(param.name, param.value, "_")
                        elif param == self.key[0]:
                            simfolder += gen_tag(param.name, param.values[i], "_") \
                                + "/"
                            datafolder += gen_tag(param.name, param.values[i], "_")
                        elif param == self.key[1]:
                            simfolder += gen_tag(param.name, param.values[j], "_") \
                                + "/"
                            datafolder += gen_tag(param.name, param.values[j], "_")

                        if k != len(self.ordered_params)-1:
                            datafolder += "_"
                        k += 1

                if self.sim_type == "cells":
                    sim = data_structures.SimulationCells(simfolder)
                else:
                    sim = data_structures.SimulationFilaments(simfolder)
                datafile = datafolder + ".txt"
                print datafile
                analysisdata = read_fnc(datafile)
                self.data[(vali, valj)] = Analysis(sim, analysisdata)

        return

    def read_data_single_instance(self, read_fnc):
        """ read data specific for single simulation instance
        --single plot type"""

        simfolder = self.simfolderbase
        datafolder = self.datafolderbase
        k = 0
        for param in self.ordered_params:
            if param.order == k:
                simfolder += gen_tag(param.name, param.value, "_") \
                    + "/"
                datafolder += gen_tag(param.name, param.value, "_")

                if k != len(self.ordered_params)-1:
                    datafolder += "_"
                k += 1

        if self.sim_type == "cells":
            sim = data_structures.SimulationCells(simfolder)
        else:
            sim = data_structures.SimulationFilaments(simfolder)

        datafile = datafolder + ".txt"
        analysisdata = read_fnc(datafile)
        self.data = Analysis(sim, analysisdata)

        return


    def set_plot_context(self):
        """ set the plot context
        like labels, title, etc"""

        self.title = gen_title(self.fixed_param)

        if self.plot_type == "2D":
            self.legend = gen_label(self.legend_param)

        elif self.plot_type == "1D":
            self.legend = gen_label(self.legend_param)

        elif self.plot_type == "phase":
            self.legend = None

        elif self.plot_type == "single":
            self.legend = None

        else:
            raise ValueError("Plot type is chosen wrongly")

        return

#############################################################################

class Analysis:

    def __init__(self, sim, data):

        self.sim = sim
        self.data = data

        return

#############################################################################

#class Separator2:
#    """ data separator for plotting purposes"""
#
#    def __init__(self, args, read_fnc):
#
#        self.sim_type = args.sim_type
#        if self.sim_type == "highdensfilaments":
#            self.legend_param, self.fixed_param = self.read_highdens_fil_param_data(args)
#        elif self.sim_type == "bidispersefilaments":
#            self.legend_param, self.fixed_param = self.read_bidisperse_fil_param_data(args)
#        elif self.sim_type == "cells":
#            self.legend_param, self.fixed_param = self.read_cell_param_data(args)
#        else:
#            raise ValueError("The chosen sim_type does not exist!")
#
#        self.out = self.read_analysis_data(args, read_fnc)
#
#        return
#
#    ##############################################################################
#
#    def read_cell_param_data(self, args):
#        """ read cell parameter data"""
#
#        legend_param = {}
#        fixed_params = {}
#        if args.kappa == -1:
#            legend_param["kappa"] = [1.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0]
#            fixed_params["eps"] = args.eps
#            fixed_params["areak"] = args.areak
#            fixed_params["fp"] = args.fp
#        elif args.fp == -1:
#            legend_param["fp"] = [0.5, 1.0, 3.0, 5.0]
#            fixed_params["eps"] = args.eps
#            fixed_params["areak"] = args.areak
#            fixed_params["kappa"] = args.kappa
#        elif args.eps == -1:
#            legend_param["eps"] = [0.05, 0.5, 1.0, 3.0, 5.0, 10.0, 20.0]
#            fixed_params["kappa"] = args.kappa
#            fixed_params["areak"] = args.areak
#            fixed_params["fp"] = args.fp
#        elif args.areak == -1:
#            legend_param["areak"] = [1.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0]
#            fixed_params["eps"] = args.eps
#            fixed_params["kappa"] = args.kappa
#            fixed_params["fp"] = args.fp
#
#        return legend_param, fixed_params
#
#    ##############################################################################
#
#    def read_highdens_fil_param_data(self, args):
#        """ read filament parameter data,
#        specific for highdens filaments data"""
#
#        legend_param = {}
#        fixed_params = {}
#        if args.kappa == -1:
#            legend_param["kappa"] = [6.0, 12.0, 24.0, 60.0]
#            fixed_params["fp"] = args.fp
#            fixed_params["density"] = args.density
#        elif args.fp == -1:
#            legend_param["fp"] = [0.001, 0.01, 0.027, 0.1]
#            fixed_params["density"] = args.density
#            fixed_params["kappa"] = args.kappa
#        elif args.density == -1:
#            legend_param["density"] = [0.8]
#            fixed_params["kappa"] = args.kappa
#            fixed_params["fp"] = args.fp
#
#        return legend_param, fixed_params
#
#    ##############################################################################
#
#    def read_bidisperse_fil_param_data(self, args):
#        """ read filament parameter data,
#        specific for bidisperse filaments data"""
#
#        legend_param = {}
#        fixed_params = {}
#        if args.kappa == -1:
#            legend_param["kappa"] = [2.5, 25.0, 400.0]
#            fixed_params["fp"] = args.fp
#            fixed_params["density"] = args.density
#        elif args.fp == -1:
#            legend_param["fp"] = [0.024, 0.24, 2.4]
#            fixed_params["density"] = args.density
#            fixed_params["kappa"] = args.kappa
#        elif args.density == -1:
#            legend_param["density"] = [0.8]
#            fixed_params["kappa"] = args.kappa
#            fixed_params["fp"] = args.fp
#
#        return legend_param, fixed_params
#
#    ##############################################################################
#
#    def gen_cell_address(self, eps, fp, areak, kappa, \
#                         analysis, dbase, analysisdbase):
#        """ data and analysis folders address generator,
#        specific for cell work"""
#
#        path1 = 'eps_' + str(eps) + \
#            '/fp_' + str(fp) + '/areak_' + str(areak)  + \
#            '/kappa_' + str(kappa)
#        path2 = analysis + '/' + analysis + \
#            '_eps_' + str(eps) + '_fp_' + str(fp) + \
#            '_areak_' + str(areak) + '_kappa_' + str(kappa) + '.txt'
#        datafolder = dbase + path1 + '/'
#        analysisfile = analysisdbase + path2
#
#        return datafolder, analysisfile
#
#    ##############################################################################
#
#    def gen_fil_address(self, density, kappa, fp, \
#                        analysis, dbase, analysisdbase):
#        """ data and analysis folders address generator,
#        specific for filament work"""
#
#        path1 = 'density_' + str(density) + '/kappa_' + str(kappa) + \
#            '/fp_' + str(fp)
#        path2 = analysis + '/' + analysis + '_density_' + str(density) + \
#            '_kappa_' + str(kappa) + '_fp_' + str(fp) + '.txt'
#        datafolder = dbase + path1 + '/'
#        analysisfile = analysisdbase + path2
#
#        return datafolder, analysisfile
#
#    ##############################################################################
#
#    def read_analysis_data(self, args, read_analysis_func):
#        """ read all the analysis data and
#        rearrange the data to make it ready for plotting,
#        namely encapsulate the data in a dictionary
#        based on the chosen parameters"""
#
#        data = {}       # carries the data per parameter set
#        sims = {}       # carries the simulation information per parameter set
#        xp = {}         # carries the x axis of the data per parameter set
#
#        param_choice = self.legend_param.keys()[0]
#
#        for p in self.legend_param[param_choice]:
#
#            if self.sim_type != "cells":
#                if param_choice == 'density':
#                    datafolder, analysisfile = self.gen_fil_address(p, \
#                        args.kappa, args.fp, args.analysisname, \
#                        args.folder, args.analysisbase)
#                elif param_choice == 'kappa':
#                    datafolder, analysisfile = self.gen_fil_address(args.density, \
#                        p, args.fp, args.analysisname, \
#                        args.folder, args.analysisbase)
#                elif param_choice == 'fp':
#                    datafolder, analysisfile = self.gen_fil_address(args.density, \
#                        args.kappa, p, args.analysisname, \
#                        args.folder, args.analysisbase)
#                sims[p] = data_structures.SimulationFilaments(datafolder, args)
#            else:
#                if param_choice == 'eps':
#                    datafolder, analysisfile = self.gen_cell_address(p, \
#                        args.fp, args.areak, args.kappa, args.analysisname, \
#                        args.folder, args.analysisbase)
#                elif param_choice == 'kappa':
#                    datafolder, analysisfile = self.gen_cell_address(args.eps, \
#                        args.fp, args.areak, p, args.analysisname, \
#                        args.folder, args.analysisbase)
#                elif param_choice == 'fp':
#                    datafolder, analysisfile = self.gen_cell_address(args.eps, \
#                        p, args.areak, args.kappa,  args.analysisname, \
#                        args.folder, args.analysisbase)
#                elif param_choice == 'areak':
#                    datafolder, analysisfile = self.gen_cell_address(args.eps, \
#                        args.fp, p, args.kappa, args.analysisname, \
#                        args.folder, args.analysisbase)
#                sims[p] = data_structures.SimulationCells(datafolder, args=args)
#
#            out = read_analysis_func(analysisfile)
#            if type(out) == tuple:
#                if len(out) == 2:
#                    x, y = out
#                    xp[p] = x
#            else:
#                y = out
#            data[p] = y
#
#        if type(out) == tuple:
#            return xp, data, param_choice, sims
#        else:
#            return self.legend_param[param_choice], data, param_choice, sims
#
#    ##############################################################################
#
#    def gen_labels_for_cells(self, param_choice, sim):
#        """ generate labels according to the parameter choice for helping in plots,
#        specific for cell work"""
#
#        eps_sym = r'$\epsilon$'
#        f_sym = r'$Pe$'
#        ka_sym = r'$\kappa_A$'
#        kb_sym = r'$\kappa_B$'
#        eps_name = "epsilon"
#        f_name = "fp"
#        ka_name = "areak"
#        kb_name = "kappa"
#
#        if param_choice == 'eps':
#            xlab = eps_sym
#            label = f_sym + "=" + "{0:.2g}".format(sim.pe) + "," + \
#                ka_sym + "=" + str(sim.areak) + "," + \
#                kb_sym + "=" + str(sim.kappa)
#            path = f_name + "=" + str(sim.fp) + "_" + \
#                ka_name + "=" + str(sim.areak) + "_" + \
#                kb_name + "=" + str(sim.kappa)
#        elif param_choice == 'fp':
#            xlab = f_sym
#            label = eps_sym + "=" + str(sim.eps) + "," + \
#                ka_sym + "=" + str(sim.areak) + "," + \
#                kb_sym + "=" + str(sim.kappa)
#            path = eps_name + "=" + str(sim.eps) + "_" + \
#                ka_name + "=" + str(sim.areak) + "_" + \
#                kb_name + "=" + str(sim.kappa)
#        elif param_choice == 'areak':
#            xlab = ka_sym
#            label = eps_sym + "=" + str(sim.eps) + "," + \
#                f_sym + "=" + "{0:.2g}".format(sim.pe) + "," + \
#                kb_sym + "=" + str(sim.kappa)
#            path = eps_name + "=" + str(sim.eps) + "_" + \
#                f_name + "=" + str(sim.fp) + "_" + \
#                kb_name + "=" + str(sim.kappa)
#        elif param_choice == 'kappa':
#            xlab = kb_sym
#            label = eps_sym + "=" + str(sim.eps) + "," + \
#                f_sym + "=" +"{0:.2g}".format(sim.pe)  + "," + \
#                ka_sym + "=" + str(sim.areak)
#            path = eps_name + "=" + str(sim.eps) + "_" + \
#                f_name + "=" + str(sim.fp) + "_" + \
#                ka_name + "=" + str(sim.areak)
#        else:
#            raise ValueError("Parameter choice is non-existing!")
#
#        return xlab, label, path
#
#    ##############################################################################
#
#    def gen_labels_for_fils(self, param_choice, sim):
#        """ generate labels according to the parameter choice for helping in plots,
#        specific for filaments works"""
#
#        dens_sym = r'$\phi$'
#        f_sym = r'$Pe$'
#        k_sym = r'$\xi/L$'
#        dens_name = "density"
#        f_name = "fp"
#        k_name = "kappa"
#
#        if param_choice == 'density':
#            xlab = dens_sym
#            label = k_sym + "=" + str(sim.xil) + "," + \
#               f_sym + "=" + "{0:.2g}".format(sim.pe)
#            path = k_name + "=" + str(sim.kappa) + "_" + \
#                f_name + "=" + str(sim.fp)
#        elif param_choice == 'fp':
#            xlab = f_sym
#            label = dens_sym + "=" + str(sim.density) + "," + \
#                k_sym + "=" + str(sim.xil)
#            path = dens_name + "=" + str(sim.density) + "_" + \
#                k_name + "=" + str(sim.kappa)
#        elif param_choice == 'kappa':
#            xlab = k_sym
#            label = dens_sym + "=" + str(sim.density) + "," + \
#                f_sym + "=" + "{0:.2g}".format(sim.pe)
#            path = dens_name + "=" + str(sim.density) + "_" + \
#                f_name + "=" + str(sim.fp)
#        else:
#            raise ValueError("Parameter choice is non-existing!")
#
#        return xlab, label, path
#
#    ##############################################################################
#
