__author__ = "Andra Stroe"
__version__ = "0.1"

import glob

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy import units as u

import constants as c


def write_emline_fits(spectrum, line):
    """ 
    Writes out the results of all the line fits for a particular source using 
    astropy
    Input:

    Return: 
        Fits table
    """
    t = Table()
    # Vacuum/laboratory properties of the line
    t['line'] = Column([line['line']], dtype='U',
                       description='Line name, which includes modifier for single or double Gaussian')
    t['wl_vacuum'] = Column([line['wl_vacuum']], unit=u.Angstrom, dtype='f', 
                            description='Wavelength in vacuum')
    t['type'] = Column([line['type']], dtype='d',
                       description='Single or double Gaussian')
    t['separation'] = Column(line['separation'], unit=u.Angstrom, dtype='f', 
                                    description='Separation from the main line, if a doublet')
    t['latex'] = Column([line['latex']], dtype='U', 
                            description='Line name written in latex format, useful for plotting purposes')
   
    # Source properties
    t['z'] = Column([spectrum.z], dtype='f', description='Source redshift, from specpro')
    
    # Redshift calculated from the line
    t['zline'] = Column([spectrum.z_line.value.value], dtype='f', 
                     unit=spectrum.z_line.value.unit, 
                     description='Redshift calculated from this emission line')
    t['zline_err'] = Column([spectrum.z_line.error.value], dtype='f', 
                     unit=spectrum.z_line.error.unit, 
                     description='Error on the redshift calculated from this emission line')
    # Redshift offset calculated between the line redshift and the redshift from 
    # specpro
    t['zoffset'] = Column([spectrum.z_offset.value.value], dtype='f', 
                     unit=spectrum.z_offset.value.unit, 
                     description='Redshift offset between line and specpro redshift')
    t['zoffset_err'] = Column([spectrum.z_offset.error.value], dtype='f', 
                     unit=spectrum.z_offset.error.unit, 
                     description='Error on the redshift offset')

    # Line fits
    # Continuum around line 
    t['cont'] =  Column([spectrum.continuum.value.value], dtype='f', 
                         unit=spectrum.continuum.value.unit, 
                         description='Continuum around line, Gaussian fit')
    t['cont_err'] =  Column([spectrum.continuum.error.value], dtype='f', 
                         unit=spectrum.continuum.error.unit, 
                         description='Error on continuum around line, Gaussian fit')
    # Gaussian fit       
    t['wl'] = Column([spectrum.wavelength.value.value], dtype='f', 
                     unit=spectrum.wavelength.value.unit, 
                     description='Restframe wavelength, Gaussian fit')
    t['wl_err'] = Column([spectrum.wavelength.error.value], dtype='f', 
                     unit=spectrum.wavelength.error.unit, 
                     description='Error on restframe wavelength, Gaussian fit')
    t['height'] = Column([spectrum.height.value.value], dtype='f', 
                     unit=spectrum.height.value.unit, 
                     description='Height, Gaussian fit')
    t['height_err'] = Column([spectrum.height.error.value], dtype='f', 
                     unit=spectrum.height.error.unit, 
                     description='Error on height, Gaussian fit')
    t['sigma'] = Column([spectrum.sigma.value.value], dtype='f', 
                     unit=spectrum.sigma.value.unit, 
                     description='Sigma, Gaussian fit')
    t['sigma_err'] = Column([spectrum.sigma.error.value], dtype='f', 
                     unit=spectrum.sigma.error.unit, 
                     description='Error on sigma, Gaussian fit')           
    t['amplitude'] = Column([spectrum.amplitude.value.value], dtype='f', 
                     unit=spectrum.amplitude.value.unit, 
                     description='Amplitude, Gaussian fit')
    t['amplitude_err'] = Column([spectrum.amplitude.error.value], dtype='f', 
                     unit=spectrum.amplitude.error.unit, 
                     description='Error on amplitude, Gaussian fit')           
    # Flux
    ff = 10.**-17
    t['flux'] = Column([spectrum.flux.value.value/ff], dtype='f', 
                     unit=spectrum.flux.value.unit*ff, 
                     description='Line flux, Gaussian fit')
    t['flux_err'] = Column([spectrum.flux.error.value/ff], dtype='f', 
                     unit=spectrum.flux.error.unit*ff, 
                     description='Error on line flux, Gaussian fit') 
    # Luminosity
    fl = 10.**40
    t['luminosity'] = Column([spectrum.luminosity.value.value/fl], dtype='f', 
                     unit=spectrum.luminosity.value.unit*fl, 
                     description='Line luminosity, Gaussian fit')
    t['luminosity_err'] = Column([spectrum.luminosity.error.value/fl], dtype='f', 
                     unit=spectrum.luminosity.error.unit*fl, 
                     description='Error on line luminosity, Gaussian fit') 
    # Equivalent width
    t['EWrest'] = Column([spectrum.ew_rest.value.value], dtype='f', 
                     unit=spectrum.ew_rest.value.unit, 
                     description='Restframe line equivalent width, Gaussian fit')
    t['EWrest_err'] = Column([spectrum.ew_rest.error.value], dtype='f', 
                     unit=spectrum.ew_rest.error.unit, 
                     description='Error on the restframe equivalent width, Gaussian fit')

    # FWHM of the line
    t['FWHM'] = Column([spectrum.fwhm.value.value], dtype='f', 
                     unit=spectrum.fwhm.value.unit, 
                     description='Restframe FWHM, not deconvolved')
    t['FWHM_err'] = Column([spectrum.fwhm.error.value], dtype='f', 
                     unit=spectrum.fwhm.error.unit, 
                     description='Error on the restframe FWHM, not deconvolved')
    # Deconvolved velocity fwhm, calculated if possible
    t['v'] = Column([spectrum.velocity_fwhm.value.value], dtype='f', 
                     unit=spectrum.velocity_fwhm.value.unit, 
                     description='Restframe, deconvolved velocity FWHM; 0 means line is unresolved')
    t['v_err'] = Column([spectrum.velocity_fwhm.error.value], dtype='f', 
                     unit=spectrum.velocity_fwhm.error.unit, 
                     description='Error on the restframe, decolvelved velocity FWHM; 99 means line is unresolved')    

    return t
    