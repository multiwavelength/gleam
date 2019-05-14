r"""
Contains a list of constants and user defines units

"""
__author__ = "Andra Stroe"
__version__ = "0.1"

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM


SN_limit = 2

# Custom units
# Flux unit
fluxu = u.erg / (u.cm ** 2 * u.s * u.Angstrom)
fluxval = 10**-16. * fluxu
fluxunit = u.def_unit('1E16 erg cm^-2 s^-1 AA^-1', fluxval, 
                        format={'latex': r'10^{16}\,erg\,cm^{-2}\,s^{-1}\,\AA^{-1}'})
ff = 10**-17 # fudge factor for flux units
fl = 10**40 # fudge factor for luminosity units

# Cosmology
H0=70
Om0=0.3
Tcmb0=2.725
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Tcmb0=Tcmb0)

# Expected line width for undetected lines in units of pixels
spectral_resolution = 5. # pixels

# Formula for resolution in the observed spectrum, not restframe!!!
# [a, b], y = a*x+b obtained from MOS_SPECTRAL_RESOLUTION.fits
m_res = 0.07445 /u.Angstrom 
c_res = 52.1895

# Atmospheric bands in their restframe
Aband = [7586., 7708.]*u.Angstrom
Bband = [6864., 6945.]*u.Angstrom
Cband = [8169., 8245.]*u.Angstrom # for Toothbrush WHT data

# Tolerance for considering emission lines part of a group and thus fitting them
# together
tolerance = 13.*u.Angstrom # Angstrom

# Range to be probed by lmfit when searching around Gaussian center hint
w = 3.*u.Angstrom #Angstrom
 
# Wavelength width used to mask lines
line_width = 20.*u.Angstrom # Angstrom

# Range used for selecting continuum left and right of the source
cont_width = 70.*u.Angstrom#*u.Angstrom
cont_plot_width = 70*u.Angstrom

# Parameters to constrain the width of the lines
pixel = 1.21*u.Angstrom # 1 pixel is 1.21 Angtrom
fwhm_min = 3*pixel
fwhm_max = 10*pixel

# Plotting parameters
# overview plot size
overview_x = 80.
overview_y = 30.

# Font sizes
# X & Y label sizes
labelsize = 20