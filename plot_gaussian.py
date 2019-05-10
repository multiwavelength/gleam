__author__ = "Andra Stroe"
__version__ = "0.1"

import os, sys
from contextlib import contextmanager

import click
import numpy as np
import matplotlibparams
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
import astropy
from astropy.table import Table

import read_files as rf
import gaussian_fitting as gf
import constants as c
import spectra_operations as so


def overplot_lines(ax, linelabels, lineloc):
    """
    Overplots emission lines on an already existing plot with axes
    Input:
        ax: matplolib axis
        linelabels: list of latex formatted names for the lines
        lineloc: list of wavelengths of the lines (most likely in Ang)
    Output:
        ax2: returns the new axis in case we want to overplot some more stuff on
             top
    """
    ax2 = ax.twiny()
    
    for xc in lineloc:
        ax.axvline(x=xc, color='0.8', linestyle='--')
 
    ax.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
    ax.patch.set_visible(False) # hide the 'canvas'
    ax2.patch.set_visible(True)
    ax2.set_xlim(ax.get_xlim())
    ax2.axes.get_xaxis().set_ticks(lineloc)
    ax2.axes.xaxis.set_ticklabels(linelabels)
    ax2.xaxis.set_label_position('top') 
    ax2.set_xticklabels(ax2.xaxis.get_majorticklabels(), rotation=90)
    return ax2



def plot_spectrum(data_path, target, spectrum, spectrum_line, spectrum_fit, 
                  line_names, line_latex, line_wls, fitted_wl, inspect, 
                  plot_fit=True, d_wl=c.cont_width):
    """
    Plots spectrum and a fit to a specific line. Overplots the entire line 
    catalogue. 
    Input:
        data_path: folder location of where the plot will be saved
        target: target structure containin info on the target
        spectrum: entire observed spectrum, format from rf.read_spectrum
        spectrum_line: observed spectrum selected around fitted lines
        spectrum_fit: fit to the spectrum (sum of continuum and Gaussians)
        list_names: list of names of the lines that were fit
        line_latex: list of latex formatted line names to overplot
        line_wls: list of the wavelengths of the overplotted lines
        fitted_wl: wavelength(s) of the fitted line
        inspect: do you want to inspect the data dynamically? then it plots the 
                 images with show(); otherwise it prints out the figure into a
                 png file
        plot_fit: do you want to plot only the observed spectrum or overplot the
                  fit to the spectrum
        d_wl: wavelength to the left and right of the line that we will plot
    Output:
        Figure in show() or a saved figure in an external png file
    """
    # Set the title to the plot
    title = (f"{target['Cluster']} {target['SourceNumber']} "
                 f"z={target['Redshift']:.3}")

    # Set the basename name of the png outfile
    basename = rf.naming_convention(data_path, target["Cluster"],
                                        target["SourceNumber"], "linefits", 
                                        target["Mode"])
    for i in range(len(line_names)):
        basename = f"{basename}.{line_names[i].strip()}"

    # Set up figure and axes
    fig = plt.figure(111)
    ax = fig.add_subplot(111)

    # Overplot the observed spectrum
    ax.plot(spectrum['wl_rest'], spectrum['flux'], color='gray')
    # Overplot the observed spectrum around the emission line in a bold colour
    ax.plot(spectrum_line['wl_rest'], spectrum_line['flux'], color='black')

    # plot zoom in on the fitted line
    sub_axes = inset_axes(ax, width="100%", height="100%", loc='upper left',
                          bbox_to_anchor=(0.5,1-0.35,.3,.3), 
                          bbox_transform=ax.transAxes)

    # Select the region around the emission line and overplot the line
    select = ( (spectrum['wl_rest'] < (np.average(fitted_wl)+d_wl) ) &  
                (spectrum['wl_rest'] > (np.average(fitted_wl)-d_wl) ) )             
    sub_axes.plot(spectrum['wl_rest'][select], spectrum['flux'][select], \
                  color='gray')

    # Mark the emission line
    for line in fitted_wl:
        sub_axes.axvline(x=line, color='k', linestyle='-')
    
    if plot_fit == True:
        plot_gaussian_fit(spectrum['wl_rest'][select], spectrum_fit, ax)
        plot_gaussian_fit(spectrum['wl_rest'][select], spectrum_fit, sub_axes)

    # Overplot the emission lines of reference
    ax2 = overplot_lines(ax, line_latex, line_wls)
    sub_axes.set_zorder(ax.get_zorder()+1) # put ax in front of ax2
    ax.patch.set_visible(False) # hide the 'canvas'
    sub_axes.patch.set_visible(True)

    # Add title to the zoomed in axis
    sub_axes.set_title(r"{}".format(title))
    ax.set_xlabel(r'$\boldsymbol{\lambda}$ (\AA)', fontsize=c.labelsize)
    ax.set_ylabel(r'$F_{\boldsymbol{\lambda}}$'+
               r'$\left(10^{-16} \mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-2}\right)$',
               fontsize=c.labelsize)
  
    # show or save the figure depending on the command line flag
    if inspect:
        fig.show()
        # While not keyboard press, i.e. while it is a mouse press go back into
        # the loop until you actually do press a key
        while not plt.waitforbuttonpress():
            pass
    if not inspect:
        fig.savefig(f"{basename}.png", format='png', bbox_inches='tight')
    plt.close(fig)



@contextmanager
def overview_plot(target, data_path, line_groups, spectrum, size_x=c.overview_x, 
                 size_y=c.overview_y, d_wl=c.cont_plot_width):
    """
    Overview plot of the spectrum of a single target, with zoom-in plots around
    the lines that were fit. The zoom-in plots have the individual Gaussian fits 
    of the lines overplotted. The lines grouped together are also plotted 
    together.
    Input:
        target: properties of the target as read in from rf.read_lof
        data_path: place where the plot will be saved
        line_groups: how many groups of lines are there, i.e. how many zoom-in 
                     plots do we want
        spectrum: spectrum of the source
        size

    """
    # Generate the title of the plot from information on the target
    title = (f"{target['Cluster']} {target['SourceNumber']} "
            f"z={target['Redshift']:.3}")
    # Basename that will be used in the savefile name
    basename = rf.naming_convention(data_path, target["Cluster"], 
                                    target["SourceNumber"], "linefits", 
                                    target["Mode"])
    Myxlabel = r'$\boldsymbol{\lambda}$ (\AA)'
    Myylabel = (r'$F_{\boldsymbol{\lambda}}$'+
               r'$\left(10^{-16} \mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-2}\right)$')
    # How many zoom-in plots will we have? As many as groups of lines, as the
    # lines grouped together will also be plotted together
    No_plots = len(line_groups)
    
    # Set figure and its size; add axis
    fig = plt.figure(figsize=(size_x, size_y))
    ax = fig.add_subplot(111)

    # Plot the base obseved spectrum to the main axis
    ax.plot(spectrum['wl_rest'], spectrum['flux'], color='k') 
    ylims = ax.get_ylim()


    j = 0

    # Function to dynamically add the fits for each emission lines as they are 
    # produced/come from a loop
    def plot_line(line, spectrum_fit):
        # Counter j to keep track in which subplot we should plot the emission
        # line
        nonlocal j

        # Select a small wavelength range around the emission line where to 
        # overplot the fit
        select = ((spectrum['wl_rest'] < (np.average(line["wl_vacuum"])+d_wl)) &  
                  (spectrum['wl_rest'] > (np.average(line["wl_vacuum"])-d_wl)) ) 

        # Create a zoomed in axis focusing on each line
        r"""
        # This preserves the ratios between the x and y axis because it is a
        # proper zoomed in axis
        axins = zoomed_inset_axes(ax, 2, loc='lower left',
                bbox_to_anchor=(0+j/(No_plots-1), 1.05),
                bbox_transform=ax.transAxes)
        """
        # Semi-random x-y ratio that does not match the main plot, but better
        # shows the lines
        axins = inset_axes(ax, size_x*0.8/No_plots, size_y/2., loc='lower left',
                bbox_to_anchor=(0+j/(No_plots-1), 1.05),
                bbox_transform=ax.transAxes)
        
        # Plot the observed spectrum in the zoomed-in axis
        axins.plot(spectrum['wl_rest'][select], spectrum['flux'][select], 
                   color='k')

        # Overplot the gaussian fit to the line in the zoomed-in axis
        plot_gaussian_fit(spectrum['wl_rest'][select], spectrum_fit, axins) 
        axins.set_ylim(ylims) 

        # Mark the emission lines fitted
        for l in line:
            axins.axvline(x=l['wl_vacuum'], color='gray', linestyle='-')
            axins.text(l['wl_vacuum'], axins.get_ylim()[1], f"{l['latex']}", 
                       ha='center', va='bottom')

        # if it's the first plot, i.e. the left-most plot, label the y axis
        if j==0: 
            axins.set_ylabel(Myylabel, fontsize=c.labelsize)

        mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")
        j = j+1  

    # Return the individual lines plots to the overview plot
    yield plot_line

    # Mark atmosphere absorption
    Aband = so.restframe_wl(c.Aband, target['Redshift'])
    Bband = so.restframe_wl(c.Bband, target['Redshift'])

    for band in [Aband, Bband]:
        ax.fill_between([band[0].value, band[1].value], ax.get_ylim()[0], 
                    ax.get_ylim()[1], facecolor="gray", alpha=0.3)
                
    ax.set_ylim(ylims) 

    # Set the proper labels
    ax.set_xlabel(Myxlabel, fontsize=c.labelsize)
    ax.set_ylabel(Myylabel, fontsize=c.labelsize)
    # Save the figure under the above-decided name
    plt.savefig(f"{basename}.png", format='png', bbox_inches='tight')
    plt.close()


def plot_gaussian_fit(wl, spectrum_fit, ax):
    """
    Plot the a line fit as a continuum + a Gaussian, whenever the line was 
    detected and not just an upper limit
    """
    for line_fit in spectrum_fit.lines:
        if isinstance(line_fit, gf.Line):
            gauss_part = gf.gauss_function(wl, line_fit.height.value, 
                            line_fit.wavelength.value, 
                            line_fit.sigma.value)
            fit_flux = gauss_part + spectrum_fit.continuum.value
            ax.plot(wl, fit_flux) 
