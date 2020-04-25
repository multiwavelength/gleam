r"""
Contains a list of constants and user defined units

"""
__author__ = "Andra Stroe"
__version__ = "0.1"

import yaml
from typing import Dict

from pydantic.dataclasses import dataclass
from astropy import units as u
from astropy.table import QTable
from astropy.cosmology import FlatLambdaCDM


class Quantity(u.SpecificTypeQuantity):
    """
    Validation of the types of unit for each parameter, to ensure the right type
    is being given.
    """

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        return cls(v)


class Length(Quantity):
    _equivalent_unit = u.m


class Frequency(Quantity):
    _equivalent_unit = (1 / u.s).unit


class Temperature(Quantity):
    _equivalent_unit = u.K


@dataclass
class Cosmology:
    """
    Cosmological parameters. Contains sensible default values and checks for 
    correct types.
    """

    H0: Frequency = 70 * u.km / (u.Mpc * u.s)
    Om0: float = 0.3
    Tcmb0: Temperature = 2.725 * u.K

    @property
    def cosmo(self):
        # Cosmology using Astropy method
        return FlatLambdaCDM(H0=self.H0, Om0=self.Om0, Tcmb0=self.Tcmb0)


@dataclass
class Setup:
    resolution: Length
    line_table: str
    mask_sky: bool


@dataclass
class FittingParameters:
    """
    Paramters needed for the fitting of the lines. Contains sensible default 
    values and checks for correct types.
    """

    SN_limit: float = 2
    spectral_resolution: int = 5
    tolerance: Length = 13.0 * u.Angstrom
    w: Length = 3 * u.Angstrom
    mask_width: Length = 20 * u.Angstrom
    cont_width: Length = 70 * u.Angstrom
    fwhm_min: int = 2
    fwhm_max: int = 15


@dataclass
class Constants:
    """
    Constants pertaining to cosmology, to fitting procedures and other telescope
    specific parameters. Contains sensible defaults and check for unit type
    correctness.
    """
    
    sky: str
    setups: Dict[str, Setup]
    fitting: FittingParameters = FittingParameters()
    cosmology: Cosmology = Cosmology()


def read_config(config_file) -> Constants:
    """
    Read YAML configuration file into a class. Not all parameters have to be 
    set. It not set, a parameter will be set to the default value. The class has 
    defaults that the config file will override. The unit types will also be
    checked for correctness.
    Input:
        config_file: path to YAML parameter file
    Output:
        return a Constants dataclass instance with the defaults and config 
        overrides.
    """
    config = yaml.safe_load(open(config_file).read())
    return Constants(**config)


a = read_config("gleamconfig.yaml")
