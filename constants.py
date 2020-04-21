r"""
Contains a list of constants and user defined units

"""
__author__ = "Andra Stroe"
__version__ = "0.1"


import configparser
from dataclasses import dataclass

from astropy import units as u
from astropy.table import QTable
from astropy.cosmology import FlatLambdaCDM

@dataclass
class Constants:
    SN_limit: u.Quantity = 2 
    H0: u.Quantity = 70 * u.km / (u.Mpc*u.s)
    Om0: u.Quantity = 0.3
    Tcmb0: u.Quantity = 2.725* u.K
    spectral_resolution: u.Quantity = 5
    m_res: u.Quantity = 0.07445/u.Angstrom
    c_res: u.Quantity  = 52.1895
    tolerance: u.Quantity = 13.0*u.Angstrom
    w: u.Quantity = 3*u.Angstrom
    line_width: u.Quantity = 20*u.Angstrom
    cont_width: u.Quantity = 70*u.Angstrom
    fwhm_min: u.Quantity = 2.0
    fwhm_max: u.Quantity = 15.0


def read_config(config_file) -> Constants:
    comments = ("#", "%", "//", "!", "--")
    delimiters = ('=', ':')
    parser = configparser.ConfigParser(delimiters=delimiters, comment_prefixes=comments, inline_comment_prefixes=comments)
    parser.optionxform = lambda option: option

    add_sec = '[fakesection]\n'+open(config_file).read()
    parser.read_string(add_sec)

# Cosmology
H0 = 70
Om0 = 0.3
Tcmb0 = 2.725
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Tcmb0=Tcmb0)
"""
cosmo = FlatLambdaCDM(H0=a.H0, Om0=a.Om0, Tcmb0=a.Tcmb0)


# Expected line width for undetected lines in units of pixels
spectral_resolution = 5.0  # pixels

# Formula for resolution in the observed spectrum, not restframe!!!
# [a, b], y = a*x+b obtained from MOS_SPECTRAL_RESOLUTION.fits
m_res = 0.07445 / u.Angstrom
c_res = 52.1895

# Tolerance for considering emission lines part of a group and thus fitting them
# together
tolerance = 13.0 * u.Angstrom

# Range to be probed by lmfit when searching around Gaussian center hint
w = 3.0 * u.Angstrom
# Wavelength width used to mask lines
line_width = 20.0 * u.Angstrom

# Range used for selecting continuum left and right of the source
cont_width = 70.0 * u.Angstrom

# Parameters to constrain the width of the lines
fwhm_min = 2.0  # pixel
fwhm_max = 15.0  # pixel

SKY = QTable.read('line_lists/Sky_bands.fits')