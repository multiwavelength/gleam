__author__ = "Andra Stroe"
__version__ = "0.1"

import os, sys
import random
from typing import List, Union, Iterable
from dataclasses import dataclass
from typing import TypeVar, Generic
import itertools

import numpy as np
import astropy
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from astropy.units.quantity import Quantity as Qty
from colorama import Fore
from colorama import init

import plot_gaussian as pg
import spectra_operations as so
import constants as c

import lmfit
from lmfit.models import GaussianModel, ConstantModel
from lmfit.model import ModelResult


Qty = astropy.units.quantity.Quantity
@dataclass
class RandomVariable:
    value: Qty
    error: Qty

    @staticmethod
    def from_param(param):
        return RandomVariable(value=param.value, error=param.stderr)
    
    @property
    def badness(self):
        return self.error / self.value

    @property
    def significance(self):
        return None if self.error is None else abs(self.value / self.error)

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

    def as_fits_table(self, line: Table) -> Table:
        """ 
        Writes out the results of all the line fits for a single emission line 
        (either stand alone or one that was part of a doublet, triplet etc), for
        a particular source using astropy
        Input:
            self: emission line fit parameters
            line: lab properties of the fitted emission line
        Return: 
            Fits table with all the lab and measured properties of a particular
            emission line
        """
        t = Table()
        # Vacuum/laboratory properties of the line
        t = add_lab_values(line, t)
    
        # Source properties
        t['z'] = Column([self.z], dtype='f', description='Source redshift, from specpro')
        
        # Redshift calculated from the line
        t['zline'] = Column([self.z_line.value.value], dtype='f', 
                        unit=self.z_line.value.unit, 
                        description='Redshift calculated from this emission line')
        t['zline_err'] = Column([self.z_line.error.value], dtype='f', 
                        unit=self.z_line.error.unit, 
                        description='Error on the redshift calculated from this emission line')
        # Redshift offset calculated between the line redshift and the redshift from 
        # specpro
        t['zoffset'] = Column([self.z_offset.value.value], dtype='f', 
                        unit=self.z_offset.value.unit, 
                        description='Redshift offset between line and specpro redshift')
        t['zoffset_err'] = Column([self.z_offset.error.value], dtype='f', 
                        unit=self.z_offset.error.unit, 
                        description='Error on the redshift offset')

        # Line fits
        # Continuum around line 
        t['cont'] =  Column([self.continuum.value.value], dtype='f', 
                            unit=self.continuum.value.unit, 
                            description='Continuum around line, Gaussian fit')
        t['cont_err'] =  Column([self.continuum.error.value], dtype='f', 
                            unit=self.continuum.error.unit, 
                            description='Error on continuum around line, Gaussian fit')
        # Gaussian fit       
        t['wl'] = Column([self.wavelength.value.value], dtype='f', 
                        unit=self.wavelength.value.unit, 
                        description='Restframe wavelength, Gaussian fit')
        t['wl_err'] = Column([self.wavelength.error.value], dtype='f', 
                        unit=self.wavelength.error.unit, 
                        description='Error on restframe wavelength, Gaussian fit')
        t['height'] = Column([self.height.value.value], dtype='f', 
                        unit=self.height.value.unit, 
                        description='Height, Gaussian fit')
        t['height_err'] = Column([self.height.error.value], dtype='f', 
                        unit=self.height.error.unit, 
                        description='Error on height, Gaussian fit')
        t['sigma'] = Column([self.sigma.value.value], dtype='f', 
                        unit=self.sigma.value.unit, 
                        description='Sigma, Gaussian fit')
        t['sigma_err'] = Column([self.sigma.error.value], dtype='f', 
                        unit=self.sigma.error.unit, 
                        description='Error on sigma, Gaussian fit')           
        t['amplitude'] = Column([self.amplitude.value.value], dtype='f', 
                        unit=self.amplitude.value.unit, 
                        description='Amplitude, Gaussian fit')
        t['amplitude_err'] = Column([self.amplitude.error.value], dtype='f', 
                        unit=self.amplitude.error.unit, 
                        description='Error on amplitude, Gaussian fit')           
        # Flux
        t['flux'] = Column([self.flux.value.value/c.ff], dtype='f', 
                        unit=self.flux.value.unit*c.ff, 
                        description='Line flux, Gaussian fit')
        t['flux_err'] = Column([self.flux.error.value/c.ff], dtype='f', 
                        unit=self.flux.error.unit*c.ff, 
                        description='Error on line flux, Gaussian fit') 
        # Luminosity
        t['luminosity'] = Column([self.luminosity.value.value/c.fl], dtype='f', 
                        unit=self.luminosity.value.unit*c.fl, 
                        description='Line luminosity, Gaussian fit')
        t['luminosity_err'] = Column([self.luminosity.error.value/c.fl], dtype='f', 
                        unit=self.luminosity.error.unit*c.fl, 
                        description='Error on line luminosity, Gaussian fit') 
        # Equivalent width
        t['EWrest'] = Column([self.ew_rest.value.value], dtype='f', 
                        unit=self.ew_rest.value.unit, 
                        description='Restframe line equivalent width, Gaussian fit')
        t['EWrest_err'] = Column([self.ew_rest.error.value], dtype='f', 
                        unit=self.ew_rest.error.unit, 
                        description='Error on the restframe equivalent width, Gaussian fit')

        # FWHM of the line
        t['FWHM'] = Column([self.fwhm.value.value], dtype='f', 
                        unit=self.fwhm.value.unit, 
                        description='Restframe FWHM, not deconvolved')
        t['FWHM_err'] = Column([self.fwhm.error.value], dtype='f', 
                        unit=self.fwhm.error.unit, 
                        description='Error on the restframe FWHM, not deconvolved')
        # Deconvolved velocity fwhm, calculated if possible
        t['v'] = Column([self.velocity_fwhm.value.value], dtype='f', 
                        unit=self.velocity_fwhm.value.unit, 
                        description='Restframe, deconvolved velocity FWHM; 0 means line is unresolved')
        t['v_err'] = Column([self.velocity_fwhm.error.value], dtype='f', 
                        unit=self.velocity_fwhm.error.unit, 
                        description='Error on the restframe, decolvelved velocity FWHM; 99 means line is unresolved')    

        # Detection flag
        t['detected'] = Column([True], dtype='bool', 
                        description='Detected or nondetected line')

        # Coverage flag
        t["covered"] = Column(
            [True], dtype="bool", description="Does the spectrum cover the line?"
        )

        return t


@dataclass
class NonDetection:
    amplitude: Qty
    continuum: RandomVariable
    restwl: float
    z: float

    @property
    def flux(self) -> Qty:
        # Line flux
        return self.amplitude.to(c.fluxu*u.Angstrom)

    @property
    def luminosity(self) -> Qty:
        # Line luminosity
        ld = c.cosmo.luminosity_distance(self.z).to(u.cm)
        v = 4.*np.pi*ld**2*self.flux
        return v

    def as_fits_table(self, line: Table) -> Table:
        """ 
        Writes out the results for a nondetected emission line within a line fit 
        (either stand alone or one that was part of a doublet, triplet etc), for
        a particular source using astropy
        Input:
            self: emission line fit parameters
            line: lab properties of the fitted emission line
        Return: 
            Fits table with all the lab and measured properties of a particular
            emission line
        """
        t = Table()
        # Vacuum/laboratory properties of the line
        t = add_lab_values(line, t)

        t['z'] = Column([self.z], dtype='f', description='Source redshift, from specpro')
        
        # Line fits
        # Continuum around line 
        t['cont'] =  Column([self.continuum.value.value], dtype='f', 
                            unit=self.continuum.value.unit, 
                            description='Continuum around line, Gaussian fit')
        t['cont_err'] =  Column([self.continuum.error.value], dtype='f', 
                            unit=self.continuum.error.unit, 
                            description='Error on continuum around line, Gaussian fit')
        t['amplitude'] = Column([self.amplitude.value], dtype='f', 
                        unit=self.amplitude.unit, 
                        description='Amplitude, Gaussian fit')
        # Flux
        t['flux'] = Column([self.flux.value/c.ff], dtype='f', 
                        unit=self.flux.unit*c.ff, 
                        description='Line flux, Gaussian fit')
        # Luminosity
        t['luminosity'] = Column([self.luminosity.value/c.fl], dtype='f', 
                        unit=self.luminosity.unit*c.fl, 
                        description='Line luminosity, Gaussian fit')
        
        # Detection flag
        t['detected'] = Column([False], dtype='bool', 
                        description='Detected or nondetected line')

        # Coverage flag
        t["covered"] = Column(
            [True], dtype="bool", description="Does the spectrum cover the line?"
        )
        return t


@dataclass
class NoCoverage:
    continuum: RandomVariable
    restwl: float
    z: float

    def as_fits_table(self, line: Table) -> Table:
        """ 
        Writes out the results for an emission line within a multiple line fit 
        for which there is no spectral coverage
        Input:
            self: emission line fit parameters
            line: lab properties of the fitted emission line
        Return: 
            Fits table with all the lab and measured properties of a particular
            emission line
        """
        t = Table()
        # Vacuum/laboratory properties of the line
        t = add_lab_values(line, t)

        t['z'] = Column([self.z], dtype='f', description='Source redshift, from specpro')
        
        # Line fits
        # Continuum around line 
        t['cont'] =  Column([self.continuum.value.value], dtype='f', 
                            unit=self.continuum.value.unit, 
                            description='Continuum around line, Gaussian fit')
        t['cont_err'] =  Column([self.continuum.error.value], dtype='f', 
                            unit=self.continuum.error.unit, 
                            description='Error on continuum around line, Gaussian fit')
        # Detection flag
        t['detected'] = Column([False], dtype='bool', 
                        description='Detected or nondetected line')

        # Coverage flag
        t["covered"] = Column(
            [False], dtype="bool", description="Does the spectrum cover the line?"
        )

        return t


@dataclass
class Spectrum:
    lines: List[Union[Line, NonDetection, NoCoverage]]
    continuum: RandomVariable


def add_lab_values(line: Table, t: Table) -> Table:
    """
    For writing out the table with the results of the fits, add the lab-measured
    properties of the emission lines.
    Input:
        line: table with properties of the emission line, such as name, 
              wavelength, type etc.
        t: table with fit results, for one emission line per source
    Output:
        table with added lines
    """
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
    return t


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
    return upper_limit


def is_good(model: ModelResult) -> bool:
    """
    Test whether the model provided a good fit: i.e. whether all the lines are 
    detected. If the model has no error bars because of a poor fit or if the
    amplitude of the Gaussian is measured at a significance lower than the 
    indicated limit that would be flagged as a bad fit
    """
    if model.result.errorbars==False: return False
    fitparams = model.params
    return all( 
        RandomVariable.from_param(fitparams[f'g{i}_amplitude']).significance>c.SN_limit 
        for i in range(len(model.components)-1)
       )


def subsets(length):
    """
    Generate all subsets of the line set, in case we want to fit only 1, 2 or
    many subsets of the original line list
    """
    for n in range(length, -1, -1):
        yield from itertools.combinations(range(length), n)


def model_selection(redshift, x, y, ystd, wl_line, fix_center=False, 
                 constrain_center=False, verbose=False) -> Spectrum:
    """
    Start with the model with the most components (constant + as many 
    Gaussians as emission lines) and find the simplest model that well 
    described the data. In case the model does not work (i.e. 1 line or more is
    non-detected), the search removes that component and attempts to fit a model
    with fewer Gaussians until it finds a fit. For all non-detected lines, an
    upper limit is estimated.
    Input: 
        redshift: redshift of the source
        x: usually an array of wavelengths
        y: usually an array of fluxes
        ystd: errors on y, usually a stdev of the flux
        wl_line: wavelengths of the lines to be fit; could be a single line or 
                 multiple lines; usually a singlet or a doublet
        fix_center: fix the center of the Gaussian in the fit
        constrain_center: constrain the center of the Gaussian to a narrow 
                          region around the expected position of the line
        verbose: print results of the fit
    Return:
        Spectrum with continuum, detected lines and upper limits
    """
    for wl_subset_indices in subsets(len(wl_line)):
        wl_subset = wl_line[list(wl_subset_indices)]
        model = fit_model(redshift, x, y, ystd, wl_subset, fix_center, 
                          constrain_center, verbose)
        if is_good(model):
            break
        else:
            if verbose==True:
                print(Fore.BLUE+f'NonDetection when fitting set of lines {wl_subset}, try a simpler model')
    else:
        return None

    # Calculate upper limit, by subtracting best fit model and then calculating
    # the upper limit from the rms noise on the residuals
    residual = y.value - model.eval(x=x.value)
    ul = upper_limit(residual, x.value)
    
    # Print a note when the spectrum does not cover one of the lines to be fit
    if verbose==True:
        for i, wl in enumerate(wl_line):
            if (i not in wl_subset_indices) and np.sum(
                ~so.mask_line(
                    x, wl_line[i], w=1.01 * c.spectral_resolution * so.resolution(x)
                )
            ) < c.spectral_resolution:
                print(Fore.BLUE + f"No spectral coverage on line {wl_line[i]}")
            else:
                continue
    
    # Add line measurements, nondetections and lines without coverage into a 
    # Spectrum class format
    fitparams = model.params
    return Spectrum(
        continuum=RandomVariable.from_param(fitparams["c"]) * y.unit,
        lines=[
            Line(
                wavelength=RandomVariable.from_param(
                    fitparams[f"g{wl_subset_indices.index(i)}_center"]
                )
                * x.unit,
                height=RandomVariable.from_param(
                    fitparams[f"g{wl_subset_indices.index(i)}_height"]
                )
                * y.unit,
                sigma=RandomVariable.from_param(
                    fitparams[f"g{wl_subset_indices.index(i)}_sigma"]
                )
                * x.unit,
                amplitude=RandomVariable.from_param(
                    fitparams[f"g{wl_subset_indices.index(i)}_amplitude"]
                )
                * x.unit
                * y.unit,
                fwhm=RandomVariable.from_param(
                    fitparams[f"g{wl_subset_indices.index(i)}_fwhm"]
                )
                * x.unit,
                z=redshift,
                restwl=wl_line[i],
                continuum=RandomVariable.from_param(fitparams["c"]) * y.unit,
            )
            if i in wl_subset_indices
            else NoCoverage(
                z=redshift,
                restwl=wl_line[i],
                continuum=RandomVariable.from_param(fitparams["c"]) * y.unit,
            )
            if np.sum(
                ~so.mask_line(
                    x, wl_line[i], w=1.01 * c.spectral_resolution * so.resolution(x)
                )
            )
            < c.spectral_resolution
            else NonDetection(
                amplitude=ul * x.unit * y.unit,
                z=redshift,
                restwl=wl_line[i],
                continuum=RandomVariable.from_param(fitparams["c"]) * y.unit,
            )
            for i in range(len(wl_line))
        ],
    )


def fit_model(
    redshift,
    x,
    y,
    ystd,
    wl_line: Iterable[Qty],
    fix_center=False,
    constrain_center=False,
    verbose=False,
) -> ModelResult:
    """
    Fits a number of Gaussians plus a constant continuum to the given data with 
    the package `lmfit'
    !!! Assumes the spectra are restframe !!!
    Input:
        redshift: redshift of the source
        x: usually an array of wavelengths
        y: usually an array of fluxes
        ystd: errors on y, usually a stdev of the flux
        wl_line: wavelengths of the lines to be fit; could be a single line or 
                 multiple lines; usually a singlet or a doublet
        fix_center: fix the center of the Gaussian in the fit
        constrain_center: constrain the center of the Gaussian to a narrow 
                          region around the expected position of the line
        verbose: print results of the fit
    Output:
        parameter fits of the Gaussian(s) + continuum
    """
    # make a model that is a number of Gaussians + a constant:
    model = sum((GaussianModel(prefix=f'g{i}_')
                 for i in range(len(wl_line.value))),
                ConstantModel())

    # define the model parameters, i.e. the continuum, wavelength of line we are
    # fitting, rough estimates for the sigma and height of the peak

    # fix the center of the Gaussian for sources where the normal fit does not 
    # want to work
    if fix_center==True:
        for i in range(len(wl_line)):
            model.set_param_hint(f'g{i}_center', vary=False)
    
    # Constrain the center within a small range around the expected wavelength
    # of the emission line
    elif constrain_center==True:
        for i, wl in enumerate(wl_line):
            model.set_param_hint(f'g{i}_center', value=wl.value, 
                                                 min=(wl-c.w).value, 
                                                 max=(wl+c.w).value)
    
    # If no fixing or constraining is done, then constrain the center to be
    # within the entire range selected around the line that is used for the 
    # fitting
    else:
        for i, wl in enumerate(wl_line):
            model.set_param_hint(f'g{i}_center', value=wl.value, 
                                                 min=(wl-c.cont_width).value, 
                                                 max=(wl+c.cont_width).value)
    # Constain the FWHM and the sigma to values set by the used. I chose 
    # reasonable values where the minimum is 3*channel width and 10*channel 
    # width
    for i, wl in enumerate(wl_line.value):
        model.set_param_hint(f'g{i}_fwhm', value=(c.fwhm_min+c.fwhm_max).value/2., 
                             min=c.fwhm_min.value, max=c.fwhm_max.value)
    for i, wl in enumerate(wl_line.value):
        model.set_param_hint(f'g{i}_sigma',
                             value=so.fwhm_to_sigma((c.fwhm_min.value+c.fwhm_max.value)/2.),
                             min=so.fwhm_to_sigma(c.fwhm_min.value), 
                             max=so.fwhm_to_sigma(c.fwhm_max.value))
    
    # Set the continuum to the median of the selected range
    ctr = np.median(y).value
    params = model.make_params(c=ctr, **{f'g{i}_center': wl
                                       for i, wl in enumerate(wl_line.value)})    
    # perform a least squares fit with errors taken into account as i.e. 1/sigma
    try:
        result:ModelResult = model.fit(y.value, params, x=x.value, 
                                       weights=1.0/ystd.value)
        if verbose==True: print(result.fit_report())
        fitparams = result.params
    except:
        print(Fore.RED+"No model could be fit")
        sys.exit("Error!")
    return result


def fit_lines(target, spectrum, line_list, line_groups, fix_center=False, 
              constrain_center=False, verbose=False, ignore_sky_lines=False):
    """
    Head function that goes through the list of lines and fits all lines for a 
    particular spectrum
    Input:
        target: Astropy table with details about the target
        spectrum: Astropy table with spectrum 
        list_list: Astropy  list of all the lines that will be fit
        line_groups: Astropy list of lines that will be fit, connected into 
                     groups based on proximity; The grouped lines will be fit
                     together as a sum of Gaussians
        fix_center: fix the center of the Gaussian in the fit
        constrain_center: constrain the center of the Gaussian to a narrow 
                          region around the expected position of the line
        verbose: print results of the fit and warnings
        ignore_sky_lines: do not mask the regions around the expected positions
                          of the sky lines
    Output:
        Astropy table containing details of emission lines and the Gaussian fits 
        to them
    """
    for group in line_groups:
        select_group = ( (line_list['wl_vacuum'].quantity>group.beginning) & 
                         (line_list['wl_vacuum'].quantity<group.ending) )
        if  ( (group.ending    < np.amax(spectrum['wl_rest'].quantity) ) & 
              (group.beginning > np.amin(spectrum['wl_rest'].quantity) ) ):
            spectrum_fit, spectrum_line = do_gaussian(line_list[select_group], 
                    line_list[~select_group], spectrum, target, fix_center, 
                    constrain_center, verbose, ignore_sky_lines)
            yield spectrum_fit, spectrum_line, line_list[select_group]
    

def do_gaussian(selected_lines, other_lines, spectrum, target, fix_center=False,
                constrain_center=False, verbose=False, ignore_sky_lines=False):
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
        fix_center: fix the center of the Gaussian in the fit
        constrain_center: constrain the center of the Gaussian to a narrow 
                          region around the expected position of the line
        verbose: print results of the fit and warnings
        ignore_sky_lines: do not mask the regions around the expected positions
                          of the sky lines
    Output:
        spectrum_fit: parameters of the fit around the doublet
        spectrum_line: extracted spectrum around the lines    
    """
    # Mask the region around the line
    mask_line = so.select_lines(selected_lines, other_lines, spectrum, target,
                                ignore_sky_lines)
                                
    # Create a new variable that contains the spectrum around the line
    spectrum_line = spectrum[mask_line]

    # Fit the gaussian(s)
    spectrum_fit = model_selection(target['Redshift'], spectrum_line['wl_rest'].quantity, 
                                spectrum_line['flux'].quantity, spectrum_line['stdev'].quantity, 
                                selected_lines['wl_vacuum'].quantity, fix_center,
                                constrain_center, verbose)

    return spectrum_fit, spectrum_line

