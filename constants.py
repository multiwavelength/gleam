r"""
Contains a list of constants and user defined units

"""
__author__ = "Andra Stroe"
__version__ = "0.1"

import yaml
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
    config = yaml.safe_load(open(config_file).read())
    print(config)
    data = {
        key: u.Quantity(value)
        for key, value in config.items()
        }
    return Constants(**data)

a = read_config('constants.param')
print(a)

cosmo = FlatLambdaCDM(H0=a.H0, Om0=a.Om0, Tcmb0=a.Tcmb0)

SKY = QTable.read('line_lists/Sky_bands.fits')