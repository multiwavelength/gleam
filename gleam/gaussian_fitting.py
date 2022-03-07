__author__ = "Andra Stroe"
__version__ = "0.1"

import sys
from typing import List, Union
from dataclasses import dataclass
import itertools

import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.table import QTable, Column, hstack
from astropy.units.quantity import Quantity as Qty
from colorama import Fore
from astropy.cosmology import FlatLambdaCDM

from lmfit.models import GaussianModel, ConstantModel
from lmfit.model import ModelResult

import gleam.spectra_operations as so
from gleam.constants import Length


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
    cosmo: FlatLambdaCDM
    resolution: Length

    @property
    def flux(self):
        # Line flux
        return RandomVariable(value=self.amplitude.value, error=self.amplitude.error,)

    @property
    def luminosity(self):
        # Line luminosity
        ld = self.cosmo.luminosity_distance(self.z).to(10 ** 28 * u.cm)
        v = 4.0 * np.pi * ld ** 2 * self.flux.value
        e = 4.0 * np.pi * ld ** 2 * self.flux.error
        return RandomVariable(v, e)

    @property
    def ew_rest(self):
        # Equivalent width of the emission line in the restframe
        v = self.amplitude.value / self.continuum.value
        e = np.abs(v) * np.sqrt(
            (self.continuum.error / self.continuum.value) ** 2
            + (self.amplitude.error / self.amplitude.value) ** 2
        )
        return RandomVariable(value=v, error=e)

    @property
    def z_offset(self):
        # Offset in the redshift calculated from the offset of the absolute
        # wavelength of the line compare to the one measured from the restframe
        # spectrum; Could be caused by imperfect redshift estimation (in case it
        # is systematically measured in more emission lines) or could be velocity
        # offset because of out/inflows
        v = self.wavelength.value / self.restwl - 1.0
        e = self.wavelength.error / self.restwl
        return RandomVariable(value=v, error=e)

    @property
    def z_line(self):
        # Redshift calculated from the offset of the absolute wavelength of the
        # line compare to the one measured from the restframe spectrum; Could
        # be caused by imperfect redshift estimation (in case it is
        # systematically measured in more emission lines) or could be velocity
        # offset because of out/inflows
        v = self.z_offset.value + self.z
        e = self.z_offset.error
        return RandomVariable(value=v, error=e)

    @property
    def velocity_fwhm(self):
        # Deconvolved velocity fwhm
        # Observed resolution at the restframe wavelength
        if (self.fwhm.value ** 2.0 - self.resolution ** 2.0) < 0:
            return RandomVariable(value=np.nan * u.km / u.s, error=np.nan * u.km / u.s)
        else:
            eff_fwhm = np.sqrt(self.fwhm.value ** 2.0 - self.resolution ** 2.0)
            v = const.c.to("km/s") * eff_fwhm / self.restwl
            e = const.c.to("km/s") * self.fwhm.error / self.restwl
            return RandomVariable(value=v, error=e)

    def as_fits_table(self, line: QTable) -> QTable:
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
        t = QTable()
        # Vacuum/laboratory properties of the line
        t = hstack([t, line])

        # Source properties
        t["z"] = Column(
            [self.z], dtype="f", description="Source redshift, from specpro"
        )

        # Redshift calculated from the line
        t["zline"] = Column(
            [self.z_line.value.value],
            dtype="f",
            unit=self.z_line.value.unit,
            description="Redshift calculated from this emission line",
        )
        t["zline_err"] = Column(
            [self.z_line.error.value],
            dtype="f",
            unit=self.z_line.error.unit,
            description="Error on the redshift calculated from this emission line",
        )
        # Redshift offset calculated between the line redshift and the redshift from
        # specpro
        t["zoffset"] = Column(
            [self.z_offset.value.value],
            dtype="f",
            unit=self.z_offset.value.unit,
            description="Redshift offset between line and specpro redshift",
        )
        t["zoffset_err"] = Column(
            [self.z_offset.error.value],
            dtype="f",
            unit=self.z_offset.error.unit,
            description="Error on the redshift offset",
        )

        # Line fits
        # Continuum around line
        t["cont"] = Column(
            [self.continuum.value.value],
            dtype="f",
            unit=self.continuum.value.unit,
            description="Continuum around line, Gaussian fit",
        )
        t["cont_err"] = Column(
            [self.continuum.error.value],
            dtype="f",
            unit=self.continuum.error.unit,
            description="Error on continuum around line, Gaussian fit",
        )
        # Gaussian fit
        t["wl"] = Column(
            [self.wavelength.value.value],
            dtype="f",
            unit=self.wavelength.value.unit,
            description="Restframe wavelength, Gaussian fit",
        )
        t["wl_err"] = Column(
            [self.wavelength.error.value],
            dtype="f",
            unit=self.wavelength.error.unit,
            description="Error on restframe wavelength, Gaussian fit",
        )
        t["height"] = Column(
            [self.height.value.value],
            dtype="f",
            unit=self.height.value.unit,
            description="Height, Gaussian fit",
        )
        t["height_err"] = Column(
            [self.height.error.value],
            dtype="f",
            unit=self.height.error.unit,
            description="Error on height, Gaussian fit",
        )
        t["sigma"] = Column(
            [self.sigma.value.value],
            dtype="f",
            unit=self.sigma.value.unit,
            description="Sigma, Gaussian fit",
        )
        t["sigma_err"] = Column(
            [self.sigma.error.value],
            dtype="f",
            unit=self.sigma.error.unit,
            description="Error on sigma, Gaussian fit",
        )
        t["amplitude"] = Column(
            [self.amplitude.value.value],
            dtype="f",
            unit=self.amplitude.value.unit,
            description="Amplitude, Gaussian fit",
        )
        t["amplitude_err"] = Column(
            [self.amplitude.error.value],
            dtype="f",
            unit=self.amplitude.error.unit,
            description="Error on amplitude, Gaussian fit",
        )
        # Flux
        t["flux"] = Column(
            [self.flux.value.value],
            dtype="f",
            unit=self.flux.value.unit,
            description="Line flux, Gaussian fit",
        )
        t["flux_err"] = Column(
            [self.flux.error.value],
            dtype="f",
            unit=self.flux.error.unit,
            description="Error on line flux, Gaussian fit",
        )
        # Luminosity
        t["luminosity"] = Column(
            [self.luminosity.value.value],
            dtype="f",
            unit=self.luminosity.value.unit,
            description="Line luminosity, Gaussian fit",
        )
        t["luminosity_err"] = Column(
            [self.luminosity.error.value],
            dtype="f",
            unit=self.luminosity.error.unit,
            description="Error on line luminosity, Gaussian fit",
        )
        # Equivalent width
        t["EWrest"] = Column(
            [self.ew_rest.value.value],
            dtype="f",
            unit=self.ew_rest.value.unit,
            description="Restframe line equivalent width, Gaussian fit",
        )
        t["EWrest_err"] = Column(
            [self.ew_rest.error.value],
            dtype="f",
            unit=self.ew_rest.error.unit,
            description="Error on the restframe equivalent width, Gaussian fit",
        )

        # FWHM of the line
        t["FWHM"] = Column(
            [self.fwhm.value.value],
            dtype="f",
            unit=self.fwhm.value.unit,
            description="Restframe FWHM, not deconvolved",
        )
        t["FWHM_err"] = Column(
            [self.fwhm.error.value],
            dtype="f",
            unit=self.fwhm.error.unit,
            description="Error on the restframe FWHM, not deconvolved",
        )
        # Deconvolved velocity fwhm, calculated if possible
        t["v"] = Column(
            [self.velocity_fwhm.value.value],
            dtype="f",
            unit=self.velocity_fwhm.value.unit,
            description="Restframe, deconvolved velocity FWHM; 0 means line is unresolved",
        )
        t["v_err"] = Column(
            [self.velocity_fwhm.error.value],
            dtype="f",
            unit=self.velocity_fwhm.error.unit,
            description="Error on the restframe, deconvolved velocity FWHM; 99 means line is unresolved",
        )

        # Detection flag
        t["detected"] = Column(
            [True], dtype="bool", description="Detected or nondetected line"
        )

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
    cosmo: FlatLambdaCDM

    @property
    def flux(self) -> Qty:
        # Line flux
        return self.amplitude

    @property
    def luminosity(self) -> Qty:
        # Line luminosity
        ld = self.cosmo.luminosity_distance(self.z).to(10 ** 28 * u.cm)
        v = 4.0 * np.pi * ld ** 2 * self.flux
        return v

    def as_fits_table(self, line: QTable) -> QTable:
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
        t = QTable()
        # Vacuum/laboratory properties of the line
        t = hstack([t, line])

        t["z"] = Column(
            [self.z], dtype="f", description="Source redshift, from specpro"
        )

        # Line fits
        # Continuum around line
        t["cont"] = Column(
            [self.continuum.value.value],
            dtype="f",
            unit=self.continuum.value.unit,
            description="Continuum around line, Gaussian fit",
        )
        t["cont_err"] = Column(
            [self.continuum.error.value],
            dtype="f",
            unit=self.continuum.error.unit,
            description="Error on continuum around line, Gaussian fit",
        )
        t["amplitude"] = Column(
            [self.amplitude.value],
            dtype="f",
            unit=self.amplitude.unit,
            description="Amplitude, Gaussian fit",
        )
        # Flux
        t["flux"] = Column(
            [self.flux.value],
            dtype="f",
            unit=self.flux.unit,
            description="Line flux, Gaussian fit",
        )
        # Luminosity
        t["luminosity"] = Column(
            [self.luminosity.value],
            dtype="f",
            unit=self.luminosity.unit,
            description="Line luminosity, Gaussian fit",
        )

        # Detection flag
        t["detected"] = Column(
            [False], dtype="bool", description="Detected or nondetected line"
        )

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

    def as_fits_table(self, line: QTable) -> QTable:
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
        t = QTable()
        # Vacuum/laboratory properties of the line
        t = hstack([t, line])

        t["z"] = Column(
            [self.z], dtype="f", description="Source redshift, from specpro"
        )

        # Line fits
        # Continuum around line
        t["cont"] = Column(
            [self.continuum.value.value],
            dtype="f",
            unit=self.continuum.value.unit,
            description="Continuum around line, Gaussian fit",
        )
        t["cont_err"] = Column(
            [self.continuum.error.value],
            dtype="f",
            unit=self.continuum.error.unit,
            description="Error on continuum around line, Gaussian fit",
        )
        # Detection flag
        t["detected"] = Column(
            [False], dtype="bool", description="Detected or nondetected line"
        )

        # Coverage flag
        t["covered"] = Column(
            [False], dtype="bool", description="Does the spectrum cover the line?"
        )

        return t


@dataclass
class Spectrum:
    lines: List[Union[Line, NonDetection, NoCoverage]]
    continuum: RandomVariable


def gauss_function(x, h, x0, sigma):
    """
    Returns the 1D Gaussian over a given range
    Input:
        x: range over which the function will be calculated
        h: height of the Gaussian
        x0: center on x of the Gaussian
        sigma: sigma of the Gaussian
    Output:
        1D Gaussian function
    """
    return h * np.exp(-((x - x0) ** 2) / (2.0 * sigma ** 2.0))


def upper_limit(y, x, SN_limit, rest_spectral_resolution):
    """
    In the case where we do not get a SN*sigma detection, we define a SN*sigma 
    upper limit as per the formula:
        SNlimit *sigma_RMS_channel* sqrt( channel_width * line_width )
    Input:
        y: spectrum in units of flux or such
        x: wavelength
        SN_limit: signal to noise limit for detections
        rest_spectral_resolution: restframed FWHM of the instrument
    """
    # Get width of a pixel in the region of the spectrum probed; in units of
    # Angstrom or such. This is pixel width/dispersion and not the resolution
    pixel = so.dispersion(x)

    # Takes SN limit into account
    upper_limit = (
        SN_limit * so.spectrum_rms(y) * np.sqrt(pixel * rest_spectral_resolution)
    )
    return upper_limit


def is_good(model: ModelResult, SN_limit) -> bool:
    """
    Test whether the model provided a good fit: i.e. whether all the lines are 
    detected. If the model has no error bars because of a poor fit or if the
    amplitude of the Gaussian is measured at a significance lower than the 
    indicated limit that would be flagged as a bad fit
    """
    if model.result.errorbars == False:
        return False
    fitparams = model.params
    return all(
        RandomVariable.from_param(fitparams[f"g{i}_amplitude"]).significance > SN_limit
        for i in range(len(model.components) - 1)
    )


def subsets(length):
    """
    Generate all subsets of the line set, in case we want to fit only 1, 2 or
    many subsets of the original line list
    """
    for n in range(length, -1, -1):
        yield from itertools.combinations(range(length), n)


def model_selection(
    redshift,
    x,
    y,
    ystd,
    wl_line,
    center_constraint,
    verbose,
    cont_width,
    w,
    SN_limit,
    rest_spectral_resolution,
    cosmo,
) -> Spectrum:
    """
    Start with the model with the most components (constant + as many 
    Gaussians as emission/absorption lines) and find the simplest model that well 
    describes the data. In case the model does not work (i.e. 1 line or more is
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
        center_constraint: fix, constrain or let free the centers of the Gaussians
                           when fitting
        verbose: print full lmfit output
        sky: sky contaminated areas to mask. If None, nothing will be masked
        cont_width: the amount left and right of the lines used for continuum
                    estimation
        mask_width: if another line B falls within the region selected for 
                    fitting line A, how big of a wavelength region should be 
                    covered for B
        w: region to be probed left and right of the starting wavelength solution
           when fitting the Gaussian
        SN_limit: signal to noise limit for detections
        rest_spectral_resolution: restframed FWHM of the instrument
        cosmo: cosmological parameters
    Return:
        Spectrum with continuum, detected lines and upper limits
    """
    for wl_subset_indices in subsets(len(wl_line)):
        wl_subset = wl_line[list(wl_subset_indices)]
        model = fit_model(
            x,
            y,
            ystd,
            wl_subset,
            center_constraint,
            verbose,
            cont_width,
            w,
            rest_spectral_resolution,
        )
        if is_good(model, SN_limit):
            break
        else:
            if verbose == True:
                print(
                    Fore.BLUE
                    + f"NonDetection when fitting set of lines {wl_subset}, try a simpler model"
                )
    else:
        return None

    # Calculate upper limit, by subtracting best fit model and then calculating
    # the upper limit from the rms noise on the residuals
    residual = y.value - model.eval(x=x.value)
    ul = upper_limit(residual, x.value, SN_limit, rest_spectral_resolution.value)

    # Print a note when the spectrum does not cover one of the lines to be fit
    if verbose == True:
        for i, wl in enumerate(wl_line):
            if (i not in wl_subset_indices) and np.sum(
                ~so.mask_line(x, wl_line[i], 1.01 * rest_spectral_resolution,)
            ) < rest_spectral_resolution // so.dispersion(x):
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
                cosmo=cosmo,
                resolution=rest_spectral_resolution,
            )
            if i in wl_subset_indices
            else NoCoverage(
                z=redshift,
                restwl=wl_line[i],
                continuum=RandomVariable.from_param(fitparams["c"]) * y.unit,
            )
            if np.sum(~so.mask_line(x, wl_line[i], 1.01 * rest_spectral_resolution,))
            < rest_spectral_resolution // so.dispersion(x)
            else NonDetection(
                amplitude=ul * x.unit * y.unit,
                z=redshift,
                restwl=wl_line[i],
                continuum=RandomVariable.from_param(fitparams["c"]) * y.unit,
                cosmo=cosmo,
            )
            for i in range(len(wl_line))
        ],
    )


def fit_model(
    x,
    y,
    ystd,
    wl_line,
    center_constraint,
    verbose,
    cont_width,
    w,
    rest_spectral_resolution,
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
        center_constraint: fix, constrain or let free the centers of the Gaussians
                           when fitting
        verbose: print full lmfit output
        sky: sky contaminated areas to mask. If None, nothing will be masked
        cont_width: the amount left and right of the lines used for continuum
                    estimation
        mask_width: if another line B falls within the region selected for 
                    fitting line A, how big of a wavelength region should be 
                    covered for B
        w: region to be probed left and right of the starting wavelength solution
           when fitting the Gaussian
        SN_limit: signal to noise limit for detections
        rest_spectral_resolution: restframed FWHM of the instrument
        cosmo: cosmological parameters
    Output:
        parameter fits of the Gaussian(s) + continuum
    """
    x = x.astype(dtype=np.float64)
    y = y.astype(dtype=np.float64)
    ystd = ystd.astype(dtype=np.float64)
    # make a model that is a number of Gaussians + a constant:
    model = sum(
        (GaussianModel(prefix=f"g{i}_") for i in range(len(wl_line.to(x.unit).value))),
        ConstantModel(),
    )

    # rescale the flux scale to get numbers comparable to the wavelength and
    # avoid numerical instabilities
    flux_scale = 1.0 / np.std(y).value * cont_width.to(x.unit).value
    y = y * flux_scale
    ystd = ystd * flux_scale

    # define the model parameters, i.e. the continuum, wavelength of line we are
    # fitting, rough estimates for the sigma and height of the peak

    # fix the center of the Gaussian for sources where the normal fit does not
    # want to work
    if center_constraint == "fixed":
        for i in range(len(wl_line)):
            model.set_param_hint(f"g{i}_center", vary=False)

    # Constrain the center within a small range around the expected wavelength
    # of the emission line
    elif center_constraint == "constrained":
        for i, wl in enumerate(wl_line):
            model.set_param_hint(
                f"g{i}_center", value=wl.to(x.unit).value, min=(wl - w).to(x.unit).value, max=(wl + w).to(x.unit).value,
            )

    # If no fixing or constraining is done, then constrain the center to be
    # within the entire range selected around the line that is used for the
    # fitting
    else:
        for i, wl in enumerate(wl_line):
            model.set_param_hint(
                f"g{i}_center",
                value=wl.to(x.unit).value,
                min=(wl - cont_width).to(x.unit).value,
                max=(wl + cont_width).to(x.unit).value,
            )

    for i, wl in enumerate(wl_line.to(x.unit).value):
        # FWHM & sigma: average between minimum and maximum expected width
        model.set_param_hint(
            f"g{i}_fwhm",
            value=rest_spectral_resolution.to(x.unit).value,
            min=rest_spectral_resolution.to(x.unit).value / 4,
        )
        model.set_param_hint(
            f"g{i}_sigma",
            value=so.fwhm_to_sigma(rest_spectral_resolution.to(x.unit).value),
            min=so.fwhm_to_sigma(rest_spectral_resolution.to(x.unit).value / 4),
        )
        # Height & amplitude: maximum y value - median of continuum
        model.set_param_hint(
            f"g{i}_height", value=max(y.value, key=abs) - np.median(y).value
        )
        model.set_param_hint(
            f"g{i}_amplitude",
            value=so.height_to_amplitude(
                max(y.value, key=abs) - np.median(y).value,
                so.fwhm_to_sigma(rest_spectral_resolution.to(x.unit).value),
            ),
        )

    # Set the continuum to the median of the selected range
    ctr = np.median(y).value
    params = model.make_params(
        c=ctr, **{f"g{i}_center": wl for i, wl in enumerate(wl_line.to(x.unit).value)}
    )

    # perform a least squares fit with errors taken into account as i.e. 1/sigma
    try:
        result: ModelResult = model.fit(
            y.value, params, x=x.value, weights=1.0 / ystd.value
        )
        if verbose == True:
            print(result.fit_report())
        try:
            result.params["c"].value = result.params["c"].value / flux_scale
        except:
            pass
        try:
            result.params["c"].stderr = result.params["c"].stderr / flux_scale
        except:
            pass
        for i, wl in enumerate(wl_line.value):
            try:
                result.params[f"g{i}_amplitude"].value = (
                    result.params[f"g{i}_amplitude"].value / flux_scale
                )
            except:
                pass
            try:
                result.params[f"g{i}_amplitude"].stderr = (
                    result.params[f"g{i}_amplitude"].stderr / flux_scale
                )
            except:
                pass
            try:
                result.params[f"g{i}_height"].value = (
                    result.params[f"g{i}_height"].value / flux_scale
                )
            except:
                pass
            try:
                result.params[f"g{i}_height"].stderr = (
                    result.params[f"g{i}_height"].stderr / flux_scale
                )
            except:
                pass

    except:
        print(Fore.RED + "No model could be fit")
        sys.exit("Error!")
    return result


def fit_lines(
    target,
    spectrum,
    line_list,
    line_groups,
    center_constraint,
    verbose,
    sky,
    cont_width,
    mask_width,
    w,
    SN_limit,
    rest_spectral_resolution,
    cosmo,
):
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
        center_constraint: fix, constrain or let free the centers of the Gaussians
                           when fitting
        verbose: print full lmfit output
        sky: sky contaminated areas to mask. If None, nothing will be masked
        cont_width: the amount left and right of the lines used for continuum
                    estimation
        mask_width: if another line B falls within the region selected for 
                    fitting line A, how big of a wavelength region should be 
                    covered for B
        w: region to be probed left and right of the starting wavelength solution
           when fitting the Gaussian
        SN_limit: signal to noise limit for detections
        rest_spectral_resolution: restframed FWHM of the instrument
        cosmo: cosmological parameters
    Output:
        Astropy table containing details of emission lines and the Gaussian fits 
        to them
    """
    for group in line_groups:
        select_group = (line_list["wavelength"] > group.beginning) & (
            line_list["wavelength"] < group.ending
        )
        if (
            (line_list["wavelength"][select_group] < np.amin(spectrum["wl_rest"])).all()
            == True
        ) | (
            (line_list["wavelength"][select_group] > np.amax(spectrum["wl_rest"])).all()
            == True
        ):
            continue
        else:
            spectrum_fit, spectrum_line = do_gaussian(
                line_list[select_group],
                line_list[~select_group],
                spectrum,
                target,
                center_constraint,
                verbose,
                sky,
                cont_width,
                mask_width,
                w,
                SN_limit,
                rest_spectral_resolution,
                cosmo,
            )
            yield spectrum_fit, spectrum_line, line_list[select_group]


def do_gaussian(
    selected_lines,
    other_lines,
    spectrum,
    target,
    center_constraint,
    verbose,
    sky,
    cont_width,
    mask_width,
    w,
    SN_limit,
    rest_spectral_resolution,
    cosmo,
):
    """
    Selects the spectrum around an emission line of interest. Then fits a single
    Gaussian plus a constant continuum to the given data with the package `lmfit'
    !!! Assumes the spectrum is restframe !!!
    Input:
        selected_lines: lines to be fit jointly
        other_lines: other lines in the catalog. Is used to mask other lines 
                     that might contaminate the continuum estimation.
        spectrum: fits table containing, the wavelength, flux and error of flux
        target: fits row with properties of the source, esp. its redshift
        center_constraint: fix, constrain or let free the centers of the Gaussians
                           when fitting
        verbose: print full lmfit output
        sky: sky contaminated areas to mask. If None, nothing will be masked
        cont_width: the amount left and right of the lines used for continuum
                    estimation
        mask_width: if another line B falls within the region selected for 
                    fitting line A, how big of a wavelength region should be 
                    covered for B
        w: region to be probed left and right of the starting wavelength solution
           when fitting the Gaussian
        SN_limit: signal to noise limit for detections
        rest_spectral_resolution: restframed FWHM of the instrument
        cosmo: cosmological parameters
    Output:
        spectrum_fit: parameters of the fit around the doublet
        spectrum_line: extracted spectrum around the lines    
    """
    # Mask the region around the line
    mask_line = so.select_lines(
        selected_lines, other_lines, spectrum, target, sky, cont_width, mask_width,
    )

    # Create a new variable that contains the spectrum around the line
    spectrum_line = spectrum[mask_line]

    # Fit the gaussian(s)
    spectrum_fit = model_selection(
        target["Redshift"],
        spectrum_line["wl_rest"],
        spectrum_line["flux_rest"],
        spectrum_line["stdev_rest"],
        selected_lines["wavelength"],
        center_constraint,
        verbose,
        cont_width,
        w,
        SN_limit,
        rest_spectral_resolution,
        cosmo,
    )

    return spectrum_fit, spectrum_line
