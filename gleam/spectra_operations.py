__author__ = "Andra Stroe"
__version__ = "0.1"


import os, sys
from dataclasses import dataclass
from typing import List
import functools
import operator

import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.table import QTable, Column
from colorama import Fore
from colorama import init

init(autoreset=True)


def average_(x, n):
    """
    Bin an array by averaging n cells together
    Input:
        x: astropy column
        n: number of cells to average
    Return:
        average each n cells
    """
    return np.average(x.reshape((-1, n)), axis=1)


def average_err(x, n):
    """
    For binning an array by averaging n cells together, propagation of errors by
    sqrt(e1**2+e2**2+e3**2.+...+en**2.)/n
    Input:
        x: astropy column, for error
        n: number of cells to average
    Return:
        geometric mean of errors
    """
    return np.sqrt(np.average((x ** 2).reshape((-1, n)), axis=1) / n)


def sigma_to_fwhm(sigma):
    """
    Convert from Gaussian sigma to FWHM
    """
    return sigma * 2.0 * np.sqrt(2.0 * np.log(2.0))


def fwhm_to_sigma(fwhm):
    """
    Convert FWHM of 1D Gaussian to sigma
    """
    return fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def height_to_amplitude(height, sigma):
    """
    Convert height of a 1D Gaussian to the amplitude
    """
    return height * sigma * np.sqrt(2 * np.pi)


def amplitude_to_height(amplitude, sigma):
    """
    Convert amplitude of a 1D Gaussian to the height
    """
    return amplitude / (sigma * np.sqrt(2 * np.pi))


def bin_spectrum(spectrum, n=2):
    """
    Bin the spectrum by averaging n number of adjacent cells together
    Input:
        spectrum: astropy spectrum, where: x axis, usually wavelength; y axis, 
                  usually flux; yerr axis, usually stdev on flux
    Output:
        binned spectrum
    """
    tsize = len(spectrum) // n * n
    spectrum = spectrum[:tsize]
    t = QTable()

    t["wl"] = Column(
        average_(spectrum["wl"], n), unit=spectrum["wl"].unit, dtype="f"
    )  # Ang
    t["flux"] = Column(
        average_(spectrum["flux"], n), unit=spectrum["flux"].unit, dtype="f"
    )  # Ang
    t["stdev"] = Column(
        average_err(spectrum["stdev"], n), unit=spectrum["flux"].unit, dtype="f"
    )

    return t


def reject_outliers(data, m=2):
    """
    Rejects outliers at a certain number of sigma away from the mean
    Input:
        data: list of input data with possible outlies
        m: number of sigma (stdev) away at which the data should be cut
    Output: 
        data: filtered data
        mask: where the data should be masked
    """
    mask = np.abs(data - np.mean(data)) <= m * np.std(data)
    return data[mask], mask


def dispersion(wl):
    """
    Returns the dispestion (wavelength width per pixel) of a 1D spectrum
    Input:
        wl: 1D wavelength in some units (preferably Astropy QTable, such that it 
        has units attached)
    Output:
        dispersion: single minimum value for dispersion in the spectrum sampled
        
    ! Writes warning if the spectrum in non uniform
    """

    diff = wl[1:] - wl[0:-1]
    diff, mask = reject_outliers(diff, 3)
    average = np.mean(diff)
    stdev = np.std(diff)
    minimum = np.min(diff)

    if stdev / average > 10 ** -3:
        print(Fore.YELLOW + "Warning: non-constant dispersion")
    return minimum


def mask_line(wl, wl_ref, mask_width):
    """
    Masks the spectrum around the listed line, within the width specified
    Input:
        wl: spectrum to be masked; preferable has unit
        wl_ref: reference wavelength that we want to mask; preferably has unit
        mask_width: width to be used for masking the line
    Output:
        mask: mask to be applied to the spectrum such that the spectrum now has 
              the line in question masked away
    """
    wl_min = wl_ref - mask_width / 2.0
    wl_max = wl_ref + mask_width / 2.0
    mask = (wl < wl_min) | (wl > wl_max)

    return mask


def spectrum_rms(y):
    """
    Calculate the rms of a spectrum, after removing the mean value
    """
    rms = np.sqrt(np.mean((y - np.mean(y)) ** 2))
    return rms


def mask_atmosphere(wl, z, sky):
    """
    Masks the spectrum around prominent optical atmospheric absorption bands
    !!! Assumes spectrum has been restframed !!!
    Input:
        wl: rest-frame wavelength spectrum
        z: redshift of the source; used to restframe the sky lines
        sky: sky bands in QTable format
    Output:
        mask: mask to be applied to the spectrum such that the spectrum now has 
              the absorption features masked
    Note: the wavelengths of the A band and B band sky absorption areas are 
    defined in the constants file
    """
    if sky is None:
        return np.ones_like(wl.value).astype(bool)
    # mask areas of absorption due to sky
    absorption = functools.reduce(
        operator.or_,
        (
            (wl > restframe_wl(band["wavelength_min"], z))
            & (wl < restframe_wl(band["wavelength_max"], z))
            for band in sky
        ),
    )

    without_absorption = ~absorption
    return without_absorption


def restframe_wl(x, z):
    """
    Transform a given spectrum x into the restframe, given the redshift
    Input:
        x: observed wavelengths of a spectrum, in Angstrom or nm for example
        z: redshift
    Return:
        restframe spectrum
    """
    return x / (1.0 + z)


def restframe_flux(x, z):
    """
    Transform a given spectrum flux density x into the restframe, given the 
    redshift
    Input:
        x: observed flux density or standard deviation of a spectrum, in 
        erg/s/cm^2/A for example
        z: redshift
    Return:
        restframe spectrum flux density 
    """
    return x * (1.0 + z)


def add_restframe(spectrum, z):
    """
    Add column with restframe spectrum onto the 1d spectrum flux
    Input:
        spectrum: Astropy table containing the 1d spectrum of a source
        z: redshift, determined externally (e.g. specpro)
    Output:
        spectrum with the new added restframe wavelength column
    """
    spectrum.add_column(restframe_wl(spectrum["wl"], z), name="wl_rest")
    spectrum.add_column(restframe_flux(spectrum["flux"],z), name="flux_rest")
    spectrum.add_column(restframe_flux(spectrum["stdev"],z), name="stdev_rest")
    
    return spectrum


def select_singleline(wl_rest, line, cont_width):
    """
    Select the region around an emission line
    !!! Assumes the spectrum is restframed
    Input:
        wl_rest: Astropy table containing the restframe wavelength for a source,
                 preferably with unit
        line: wavelength of the line of interest, preferably with unit
        cont_width: width of the region around the lines to be selected, 
                    preferably with unit
    Output:
        mask to select region around the line of interest
    """
    wl_min = line - cont_width
    wl_max = line + cont_width
    mask = (wl_rest > wl_min) & (wl_rest < wl_max)
    return mask


def select_lines(
    selected_lines, other_lines, spectrum, target_info, sky, cont_width, mask_width,
):
    """
    Masks the spectrum in the vicinity of the line of interest. It should leave 
    unmasked the actual line and lots of continuum, but mask other neighboring 
    lines we want to fit, that might contaminate the continuum estimation
    Input: 
        selected_lines: table of all the lines to be fit
        other_lines: the other lines in the table that will be masked
        spectrum: 1d spectrum, error and wavelength, as read in by 
                  read_files.read_spectrum with extra column for the restframe 
                  wavelength
        target: ancillary information on the source, such as RA, DEC or redshift
                as produced by read_files.read_lof()
        cont_width: amount of wavelength coverage on each side of the line that
                    will be taken into account
    Output:
        wl_masked: masked wavelength coverage, with only the line of interest 
                   and the continuum; all other lines masked
    """
    z_ref = target_info["Redshift"]
    wl_rest = spectrum["wl_rest"]
    flux = spectrum["flux_rest"]
    # Restframe resolution
    res_rest = dispersion(wl_rest)
    # mask all lines, but the line we are interested in
    masked_otherlines = np.full(np.shape(wl_rest), True)
    for line in map(QTable, other_lines):
        mask = mask_line(wl_rest, line["wavelength"], mask_width)
        masked_otherlines = masked_otherlines & mask

    # select the lines of interest
    select_lines = np.full(np.shape(wl_rest), False)
    for line in map(QTable, selected_lines):
        mask = select_singleline(wl_rest, line["wavelength"], cont_width)
        select_lines = select_lines | mask

    # mask the atmospheric lines, if masking them is enabled
    masked_atm = mask_atmosphere(wl_rest, z_ref, sky)
    masked_all = masked_atm & masked_otherlines & select_lines

    return masked_all


def group_lines(line_list, tolerance):
    """
    Group together lines within a wavelength tolerance. These will be fit 
    together as a sum of Gaussian, rathen than independently.
    Input: 
        line_list: Astropy table containing lines and their properties
        t: how far from each other can lines be to still be considered a group. 
           tolerance is read from the constants file
    Output:
        Groups of lines of type Component as returned by connected components
    """
    wl = [(x - tolerance, x + tolerance) for x in line_list["wavelength"]]
    line_groups = connected_components(wl, left, right)
    return line_groups


@dataclass
class Component:
    segments: List
    beginning: float
    ending: float


# offline algorithm
def connected_components(segments, left, right):
    """
    For a set of segments (in my case wavelength segments), check whether they
    overlap and group them together.
    Input:
        segments: list of pair of the outmost edges of the segment; beginning 
                  and end of the segment
        left: function defining how to get the left-most edge
        right: function defining how to get the right-most edge
    Output:
        List of components (as defined in the Component class): grouped list of 
        each of the segments that overlap, with their overall left and right
        side edges
    """
    segments = sorted(segments, key=left)
    try:
        head, *segments = segments
    except:
        return []

    ref_component = Component(segments=[head], beginning=left(head), ending=right(head))
    components = [ref_component]

    for segment in segments:
        opening = left(segment)
        closing = right(segment)
        if ref_component.ending > opening:
            ref_component.segments.append(segment)
            if closing > ref_component.ending:
                ref_component.ending = closing
        else:
            ref_component = Component(
                segments=[segment], beginning=opening, ending=closing
            )
            components.append(ref_component)
    return components


def left(x):
    return x[0]


def right(x):
    return x[1]
