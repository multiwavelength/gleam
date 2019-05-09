__author__ = "Andra Stroe"
__version__ = "0.1"

import os, sys
import random
from typing import List
from dataclasses import dataclass
from typing import TypeVar, Generic

import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column

import plot_gaussian as pg
import spectra_operations as so
import constants as c

import lmfit
from lmfit.models import GaussianModel, ConstantModel

Qty = TypeVar("Qty")
@dataclass
class RandomVariable(Generic[Qty]):
    value: Qty
    error: Qty

    @staticmethod
    def from_param(param):
        return RandomVariable(value=param.value, error=param.stderr)
    
    @property
    def badness(self):
        return self.error / self.value

    @property
    def SN(self):
        return self.value / self.error

    @property
    def significance(self):
        return abs(self.value / self.error)

    def __mul__(self, other):
        """ self * other """
        return RandomVariable(value=self.value * other, error=self.error * other)

    def __rmul__(self, other):
        """ other * self """
        return self * other

@dataclass
class Line:
    wavelength: RandomVariable
    height: RandomVariable
    sigma: RandomVariable
    amplitude: RandomVariable
    continuum: RandomVariable
    fwhm: RandomVariable
    restwl: float
    z: float

    @property
    def flux(self):
        # Line flux
        return RandomVariable(value=self.amplitude.value.to(c.fluxu*u.Angstrom), 
                              error=self.amplitude.error.to(c.fluxu*u.Angstrom))

    @property
    def luminosity(self):
        # Line luminosity
        ld = c.cosmo.luminosity_distance(self.z).to(u.cm)
        v = 4.*np.pi*ld**2*self.flux.value
        e = 4.*np.pi*ld**2*self.flux.error
        return RandomVariable(v, e)

    @property
    def ew_rest(self):
        # Equivalent width of the emission line in the restframe
        v = self.amplitude.value/self.continuum.value
        e = v * np.sqrt( (self.continuum.error/self.continuum.value)**2 + 
                         (self.amplitude.error/self.amplitude.value)**2)
        return RandomVariable(value=v, error=e)

    @property
    def z_offset(self):
        # Offset in the redshift calculated from the offset of the absolute 
        # wavelength of the line compare to the one measured from the restframe 
        # spectrum; Could be caused by imperfect redshift estimation (in case it 
        # is systematically measured in more emission lines) or could be velocity 
        # offset because of out/inflows
        v = self.wavelength.value/self.restwl-1.
        e = self.wavelength.error/self.restwl
        return RandomVariable(value=v, error=e)

    @property
    def z_line(self):
        # Redshift calculated from the offset of the absolute wavelength of the
        # line compare to the one measured fromt the restframe spectrum; Could
        # be caused by imperfect redshift estimation (in case it is 
        # systematically measured in more emission lines) or could be velocity 
        # offset because of out/inflows
        v = self.z_offset.value+self.z
        e = self.z_offset.error
        return RandomVariable(value=v, error=e)

    @property
    def velocity_fwhm(self):
        # Deconvolved velocity fwhm
        # Observed resolving power at the line wavlength
        res_power = c.m_res*self.restwl*(1+self.z)+c.c_res
        # Observed resolution at the restframe wavelength
        resolution = self.restwl/res_power
        if (self.fwhm.value**2. - resolution**2. )<0: 
            return RandomVariable(value=0*u.km/u.s, error=99*u.km/u.s)
        else: 
            eff_fwhm = np.sqrt(self.fwhm.value**2. - resolution**2. )
            v = const.c.to('km/s')*eff_fwhm/self.restwl
            e = const.c.to('km/s')*self.fwhm.error/self.restwl
            return RandomVariable(value=v,error=e)
        

    
@dataclass
class Spectrum:
    lines: List[Line]
    continuum: RandomVariable

def gauss_function(x, a, x0, sigma):
    """
    Returns the 1D Gaussian over a given range
    Input:
        x: range over which the function will be calculated
        a: amplitude of the Gaussian
        x0: center on x of the Gaussian
        sigma: sigma of the Gaussian
    Output:
        1D Gaussian function
    """
    return a*np.exp(-(x-x0)**2/(2.*sigma**2.))


def upper_limit(y,x):
    """
    In the case where we do not get a SN*sigma detection, we define a SN*sigma 
    upper limit as per the formula:
        3*sigma_RMS_channel* sqrt( channel_width * line_width )
    Input:
        y: spectrum in units of flux or such
        x: wavelength
    """
    # Get width of a pixel in the region of the spectrum probed; in units of 
    # Angstrom or such. This is pixel width and not resolution
    pixel = so.resolution(x) 
    # Resolution is assumed to be a multiple of the pixel size, e.g. 5
    # Takes SN limit into account
    upper_limit = ( c.SN_limit * so.spectrum_rms(y) * 
                    np.sqrt(pixel**2*c.spectral_resolution) )
    return upper_limit.to(x.unit*y.unit)


def fit_gaussian(redshift, x, y, ystd, wl_line, fix_center=False, 
                 constrain_center=False, ctr=0.5, peak_sigma=5., 
                 peak_height=1.) -> Spectrum:
    """
    Fits a number of Gaussians plus a constant continuum to the given data with 
    the package `lmfit'
    !!! Assumes the spectra are restframe !!!
    Input:
        x: usually an array of wavelengths
        y: usually an array of fluxes
        ystd: errors on y, usually a stdev of the flux
        wl_line: wavelengths of the lines to be fit; could be a single line or 
                 multiple lines; usually a singlet or a doublet
        c, peak_sigma, peak_height: some reasonable guesses for the constant,
            the sigma of the peak and the amplitude of the peak
    Output:
        parameter fits of the Gaussian(s) + continuum
    """
    # make a model that is a number of Gaussians + a constant:
    model = sum((GaussianModel(prefix=f'g{i}_')
                 for i in range(len(wl_line))),
                ConstantModel())

    # define the model parameters, i.e. the continuum, wavelength of line we are
    # fitting, rough estimates for the sigma and height of the peak

    # fix the center of the Gaussian for sources where the normal fit does not 
    # want to work
    
    if fix_center==True:
        for i in range(len(wl_line)):
            model.set_param_hint(f'g{i}_center', vary=False)
    if constrain_center==True:
        for i, wl in enumerate(wl_line):
            model.set_param_hint(f'g{i}_center', value=wl, min=wl-c.w, max=wl+c.w)
    params = model.make_params(c=ctr, **{f'g{i}_center': wl
                                       for i, wl in enumerate(wl_line)})    
    # perform a least squares fit with errors taken into account as i.e. 1/sigma
    try:
        result = model.fit(y, params, x=x, weights=1.0/ystd)
    except:
        return None
    
    # Return nothing if fit was not possible
    for param in result.params.values():
        if RandomVariable.from_param(param).error == None:
            return None
        if RandomVariable.from_param(param).badness > 10:
            return None 

    params = result.params
    
    return Spectrum(continuum=RandomVariable.from_param(params['c'])*y.unit,
                    lines=[Line(wavelength=RandomVariable.from_param(params[f'g{i}_center'])*x.unit,
                                height=RandomVariable.from_param(params[f'g{i}_height'])*y.unit,
                                sigma=RandomVariable.from_param(params[f'g{i}_sigma'])*x.unit,
                                amplitude=RandomVariable.from_param(params[f'g{i}_amplitude'])*x.unit*y.unit,
                                fwhm=RandomVariable.from_param(params[f'g{i}_fwhm'])*x.unit,
                                z=redshift,
                                restwl=wl_line[i]*wl_line.unit, 
                                continuum=RandomVariable.from_param(params['c'])*y.unit)
                           for i in range(len(wl_line))])



def fit_lines(target, spectrum, line_list, line_groups, fix_center=False, 
              constrain_center=False):
    """
    Head function that goes through the list of lines and fits all lines for a 
    particular spectrum
    Input:
        target: Astropy table with details about the target
        spectrum: Astropy table with spectrum 
        list_list: Astropy  list of all the lines that will be fit
        line_groups: Astropy list of lines that will be fit, connected into 
                     groups based on proximity; The grouped lines will be fit
                     together as a sum of Gaussianss
    Output:
        Astropy table containing details of emission lines and the Gaussian fits 
        to them
    """
    for group in line_groups:
        select_group = ( (line_list['wl_vacuum']>group.beginning) & 
                         (line_list['wl_vacuum']<group.ending) )
        if  ( (group.ending    < np.amax(spectrum['wl_rest']) ) & 
              (group.beginning > np.amin(spectrum['wl_rest']) ) ):
            spectrum_fit, spectrum_line = do_gaussian(line_list[select_group], 
                    line_list[~select_group], spectrum, target, fix_center, 
                    constrain_center)
            yield spectrum_fit, spectrum_line, line_list[select_group]
    

def do_gaussian(selected_lines, other_lines, spectrum, target, fix_center=False,
                constrain_center=False):
    """
    Selects the spectrum around an emission line of interest. Then fits a single
    Gaussian plus a constant continuum to the given data with the package `lmfit'
    !!! Assumes the spectrum is restframe !!!
    Input:
        x: usually a wavelength
        y: usually a flux
        ystd: errors on y, usually a stdev of the flux
        wl_line: wavelength of the line to be fit
        c, peak_sigma, peak_amplitude: some reasonable guesses for the constant,
            the sigma of the peak and the amplitude of the peak
        fix_center: fix the center of the Gaussian
    Output:
        spectrum_fit: parameters of the fit around the doublet
        spectrum_line: extracted spectrum around the lines    
    """
    # Mask the region around the line
    mask_line = so.select_lines(selected_lines, other_lines, spectrum, target)
    # Create a new variable that contains the spectrum around the line
    spectrum_line = spectrum[mask_line]

    # Fit the gaussian(s)
    spectrum_fit = fit_gaussian(target['Redshift'], spectrum_line['wl_rest'], 
                                spectrum_line['flux'], spectrum_line['stdev'], 
                                selected_lines['wl_vacuum'], fix_center,
                                constrain_center)
    
    return spectrum_fit, spectrum_line

