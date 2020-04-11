r"""
Contains a list of constants and user defined units

"""
__author__ = "Andra Stroe"
__version__ = "0.1"

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

# Signal to noise limit for a detection
SN_limit = 2

# Custom units -- change as necessary
# Flux unit for flux calibrated data
"""
fluxu = u.erg / (u.cm ** 2 * u.s * u.Angstrom)
fluxval = 10**-16. * fluxu
fluxunit = u.def_unit('1E-16 erg cm^-2 s^-1 AA^-1', fluxval, 
                        format={'latex': r'10^{-16}\,erg\,cm^{-2}\,s^{-1}\,\AA^{-1}'})
ff = 10**-17 # fudge factor for flux units
fl = 10**40 # fudge factor for luminosity units
"""
# Flux unit for MMT where data is not flux calibrated
fluxu = u.ct
fluxval = fluxu
fluxunit = u.ct
ff = 1  # fudge factor for flux units
fl = 1  # fudge factor for luminosity units

# Cosmology
H0 = 70
Om0 = 0.3
Tcmb0 = 2.725
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Tcmb0=Tcmb0)

# Expected line width for undetected lines in units of pixels
spectral_resolution = 5.0  # pixels

# Formula for resolution in the observed spectrum, not restframe!!!
# [a, b], y = a*x+b obtained from MOS_SPECTRAL_RESOLUTION.fits
m_res = 0.07445 / u.Angstrom
c_res = 52.1895

# Atmospheric bands in their restframe
Aband = [7586.0, 7658.0] * u.Angstrom
Bband = [6864.0, 6945.0] * u.Angstrom

# Tolerance for considering emission lines part of a group and thus fitting them
# together
tolerance = 13.0 * u.Angstrom

# Range to be probed by lmfit when searching around Gaussian center hint
w = 3.0 * u.Angstrom
# Wavelength width used to mask lines
line_width = 20.0 * u.Angstrom

# Range used for selecting continuum left and right of the source
cont_width = 70.0 * u.Angstrom
cont_plot_width = 70 * u.Angstrom

# Parameters to constrain the width of the lines
fwhm_min = 2.0  # pixel
fwhm_max = 15.0  # pixel

# Plotting parameters
# overview plot size
overview_x = 60.0
overview_y = 15.0

# Font sizes
# X & Y label sizes
labelsize = 30
