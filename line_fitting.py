__author__ = "Andra Stroe"
__version__ = "0.1"

import os, sys

import numpy as np
import astropy
from astropy import units as u
from astropy.table import Table, Column
from astropy.io import fits
from colorama import Fore
from colorama import init
init(autoreset=True)

import read_files as rf
import gaussian_fitting as gf
import plot_gaussian as pg
import spectra_operations as so
import write_files as wf


def run_main(data_path, target, line_list, inspect=False, fix_center=False,
             constrain_center=False, bin1=1):
    """
    For a target/galaxy, read the spectrum and perform the line fitting for each 
    line within the list of lines
    Input:
        data_path: location on disk (folder) of the target spectrum 
        target: observed specturm in specpro format 
        line_list: list of emission lines to measure; in astropy table format
        inspect: if true, show the plots; otherwise write to disk
    Output:
        fits of emission lines and plots for each fitted lines
    """
    # Read spectrum for the current source from outside file
    spectrum = rf.read_spectrum(data_path, target["Cluster"],
                                target["SourceNumber"], target["Mode"])
    if bin1>1: spectrum = so.bin_spectrum(spectrum, bin1)

    
    
    # Given its redshift, calculate restframe spectrum
    #spectrum = so.add_restframe(spectrum, target["Redshift"])
    spectrum = so.add_restframe(spectrum, target["Redshift"])

    # Find groups of nearby lines in the input table that will be fit together
    line_groups = so.group_lines(line_list)
    with pg.overview_plot(target, data_path, line_groups, spectrum) as plot_line:
        tables = []
        # Set the name to the exported plot in png format            
        for spectrum_fit, spectrum_line, lines in gf.fit_lines(
                target, spectrum, line_list, line_groups, fix_center, 
                constrain_center):
            # Make a plot/fit a spectrum if the line in within the rest-frame 
            # spectral coverage of the source
            #pg.line_centres(line_list)
            pg.plot_spectrum(data_path, target, spectrum, spectrum_line,
                             spectrum_fit, lines["line"], 
                             line_list["latex"],
                             line_list["wl_vacuum"],
                             lines["wl_vacuum"],
                             inspect)
            if spectrum_fit is not None: 
                for (line_fit, line)  in zip(spectrum_fit.lines, lines):
                    tables.append(wf.write_emline_fits(line_fit, line))                       
                plot_line(lines, spectrum_fit)

    try:
        outtable = astropy.table.vstack(tables)
        outfile = '{}.fits'.format(rf.naming_convention(data_path, 
                                target['Cluster'], target['SourceNumber'], 
                                'linefits', target["Mode"]))
        outtable.write(outfile, overwrite=True)
    except:
        print(Fore.RED+"Warning: no emission line fits in source: {} {}".format(
                        target["Cluster"], target["SourceNumber"]))

    
