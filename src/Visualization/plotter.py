#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 16:31:55 2017

@author: duman
"""

""" Data types and functions for plotting"""

###

##############################################################################

import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import read_write
import misc_tools
import data_separator

import seaborn as sns
sns.set(style="white",context='paper',
        font_scale=1.2,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 30,
            'font.family': 'sans',"figure.dpi":300,
            "xtick.major.size": 8, "ytick.major.size": 8,
            'grid.linestyle': '--'})

##############################################################################

class Subplots:
    """ subplots structure"""

    totcnt = -1             # Total number of subplots -- static member

    def __init__(self, f, l=1.0, s=0.0, b=0.0, t=1):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box
        self.sep = s        # Separation distance between subplots
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction

        return

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

class GeneralPlot(Subplots):
    """ general plot structure"""

    def __init__(self, num_ticks=5,
                 ax_len=1.0, ax_b=0.0,
                 ax_sep=0.0, total_subplots_in_x=1):

        self.fig = plt.figure()
        Subplots.__init__(self, self.fig, \
                          ax_len, ax_sep, ax_b, total_subplots_in_x)
        self.ax0 = self.addSubplot()

        return

    def plot_1d(self, x, y, sims, savebase, analysisname, param_choice, separator, \
                props, savepdf):
        """ plot 1D analysis data with 1 control parameter as legend,
        and 1 control parameter in the x axis"""

        os.system("mkdir -p " + savebase)
        savepath = savebase + analysisname + "/"
        os.system("mkdir -p " + savepath)
        savepath = savepath + param_choice + "/"
        os.system("mkdir -p " + savepath)

        keys = np.sort(sims.keys())
        xp = np.zeros((len(keys)), dtype=np.float32)
        yp = np.zeros((len(keys)), dtype=np.float32)
        for j, key in enumerate(keys):
            xp[j] = x[j]
            yp[j] = y[key]
            sim = sims[key]

            xlab, title, path = separator.gen_labels_for_cells(param_choice, sim)
            print xp[j], yp[j]
        line0 = self.ax0.scatter(xp, yp, \
                         s=100)

        ### save address

        savepath = savepath + analysisname + "_per_" + param_choice + "_" + path
        print savepath

        ### scales

        self.ax0.set_xscale(props.xscale)
        self.ax0.set_yscale(props.yscale)

        ### title

        self.ax0.set_title(title, fontsize=30)

        ### labels

        self.ax0.set_xlabel(xlab, fontsize=40)
        self.ax0.set_ylabel(props.ylab, fontsize=40)

        ### limits
        if props.set_xlim:
            self.ax0.set_xlim(props.xlim)
        if props.set_ylim:
           self.ax0.set_ylim(props.ylim)

        ### ticks

        if props.set_xticks:
            self.ax0.xaxis.set_ticks(props.xticks)
        if props.set_yticks:
            self.ax0.yaxis.set_ticks(props.yticks)
        self.ax0.tick_params(axis='both', which='major', labelsize=30)

        ### legend

        #self.ax0.legend(bbox_to_anchor=(0.005, 0.,0.65, 1.), loc=2, borderaxespad=0., \
        #    prop={'size': 20}, mode="expand", frameon=False)

        ### save

        if savepdf:
            plt.savefig(savepath+".pdf", dpi=300, bbox_inches='tight', pad_inches=0.08)
        else:
            plt.savefig(savepath+".png", dpi=300, bbox_inches='tight', pad_inches=0.08)

        self.fig.clf()
        plt.close()

        return

    def plot_2d(self, x, y, sims, savebase, analysisname, param_choice, separator, \
                props, savepdf):
        """ plot 2D analysis data with 1 control parameter as legend"""

        os.system("mkdir -p " + savebase)
        savepath = savebase + analysisname + "/"
        os.system("mkdir -p " + savepath)
        savepath = savepath + param_choice + "/"
        os.system("mkdir -p " + savepath)

        keys = np.sort(sims.keys())
        for key in keys:
            xp = x[key]
            yp = y[key]
            sim = sims[key]

            legend, title, path = separator.gen_labels_for_cells(param_choice, sim)

            label = legend + "=" + "{0:.2g}".format(key)
            line0 = self.ax0.plot(xp, yp, \
                                 linewidth=2.0, label=label)

        ### save address

        savepath = savepath + analysisname + "_per_" + param_choice + "_" + path
        print savepath

        ### scales

        self.ax0.set_xscale(props.xscale)
        self.ax0.set_yscale(props.yscale)

        ### title

        self.ax0.set_title(title, fontsize=30)

        ### labels

        self.ax0.set_xlabel(props.xlab, fontsize=40)
        self.ax0.set_ylabel(props.ylab, fontsize=40)

        ### limits

        if props.set_xlim:
            self.ax0.set_xlim(props.xlim)
        if props.set_ylim:
           self.ax0.set_ylim(props.ylim)

        ### ticks

        if props.set_xticks:
            self.ax0.xaxis.set_ticks(props.xticks)
        if props.set_yticks:
            self.ax0.yaxis.set_ticks(props.yticks)
        self.ax0.tick_params(axis='both', which='major', labelsize=30)

        ### legend

        self.ax0.legend(bbox_to_anchor=(0.005, 0.,0.65, 1.), loc=2, borderaxespad=0., \
            prop={'size': 20}, mode="expand", frameon=False)

        ### save

        if savepdf:
            plt.savefig(savepath+".pdf", dpi=300, bbox_inches='tight', pad_inches=0.08)
        else:
            plt.savefig(savepath+".png", dpi=300, bbox_inches='tight', pad_inches=0.08)

        self.fig.clf()
        plt.close()

        return

##############################################################################

def plot_2d(plot_specific):
    """ plot analysis data in 2d plot type"""

    fig = plt.figure()
    subp = Subplots(fig)
    ax0 = subp.addSubplot()

    def func_wrapper(data, sep, savepdf):

        ### generate the save folder path to store the figure

        savefolder = sep.savefolderbase + \
            sep.legend_param.physicalparamname.upper() + "/"
        os.system("mkdir -p " + savefolder)

        ### plot the figure

        updated_fig = plot_specific(data, sep, savepdf, fig, ax0)

        ### generate the save file path address to store the figure

        data_separator.assign_physicalvalues(sep.fixed_param, data[data.keys()[0]].sim)
        savefilepath = savefolder + data_separator.gen_path(sep.fixed_param)
        if savepdf:
            savefilepath += ".pdf"
        else:
            savefilepath += ".png"
        print "Saving the figure: ", savefilepath

        ### store the figure

        plt.savefig(savefilepath, dpi=300, bbox_inches='tight', pad_inches=0.08)
        fig.clf()

        return

    return func_wrapper

##############################################################################

def plot_1d(plot_specific):
    """ plot analysis data in 1d plot type"""

    fig = plt.figure()
    subp = Subplots(fig)
    ax0 = subp.addSubplot()

    def func_wrapper(data, sep, savepdf):

        ### generate the save folder path to store the figure

        savefolder = sep.savefolderbase + \
            sep.legend_param.physicalparamname.upper() + "/"
        os.system("mkdir -p " + savefolder)

        ### plot the figure

        updated_fig = plot_specific(data, sep, savepdf, fig, ax0)

        ### generate the save file path address to store the figure

        data_separator.assign_physicalvalues(sep.fixed_param, data[data.keys()[0]].sim)
        savefilepath = savefolder + \
            sep.control_param.physicalparamname.upper() + "_" + data_separator.gen_path(sep.fixed_param)
        if savepdf:
            savefilepath += ".pdf"
        else:
            savefilepath += ".png"
        print "Saving the figure: ", savefilepath

        ### store the figure

        plt.savefig(savefilepath, dpi=300, bbox_inches='tight', pad_inches=0.08)
        fig.clf()

        return

    return func_wrapper

##############################################################################

def plot_phase(plot_specific):
    """ plot analysis data in phase plot type"""

    fig = plt.figure()
    subp = Subplots(fig)
    ax0 = subp.addSubplot()

    def func_wrapper(data, sep, savepdf):

        ### generate the save folder path to store the figure

        savefolder = sep.savefolderbase + \
            sep.control_param[0].physicalparamname.upper() + "/"
        os.system("mkdir -p " + savefolder)

        ### plot the figure

        updated_fig = plot_specific(data, sep, savepdf, fig, ax0)

        ### generate the save file path address to store the figure

        data_separator.assign_physicalvalues(sep.fixed_param, data[data.keys()[0]].sim)
        savefilepath = savefolder + \
            sep.control_param[1].physicalparamname.upper() + "_" + data_separator.gen_path(sep.fixed_param)
        if savepdf:
            savefilepath += ".pdf"
        else:
            savefilepath += ".png"
        print "Saving the figure: ", savefilepath

        ### store the figure

        plt.savefig(savefilepath, dpi=300, bbox_inches='tight', pad_inches=0.08)
        fig.clf()

        return

    return func_wrapper

##############################################################################

