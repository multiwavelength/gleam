r"""
Contains a list of constants and user defined units

"""
__author__ = "Andra Stroe"
__version__ = "0.1"

import yaml
from dataclasses import replace, asdict, is_dataclass, field
from typing import Dict, Optional, List, Union, Literal
import operator
import functools


from pydantic.dataclasses import dataclass
from astropy import units as u
from astropy.table import QTable
from astropy.cosmology import FlatLambdaCDM


AllLines = Literal["all"]


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


def override(basic, overrides):
    type_b = type(basic)
    if overrides is None:
        return basic
    if is_dataclass(basic):
        basic = asdict(basic)
    if is_dataclass(overrides):
        overrides = asdict(overrides)
    return type_b(
        **{
            **basic,
            **{key: value for key, value in overrides.items() if value is not None},
            **{
                key: override(value, overrides.get(key))
                for key, value in basic.items()
                if isinstance(value, dict)
            },
        }
    )


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
class FittingParametersOverrides:
    """
    Paramters needed for the fitting of the lines. Contains sensible default 
    values and checks for correct types.
    """

    SN_limit: Optional[float] = None
    tolerance: Optional[Length] = None
    w: Optional[Length] = None
    mask_width: Optional[Length] = None
    cont_width: Optional[Length] = None
    fwhm_min: Optional[int] = None
    fwhm_max: Optional[int] = None


@dataclass
class FittingParameters:
    """
    Paramters needed for the fitting of the lines. Contains sensible default 
    values and checks for correct types.
    """

    SN_limit: float = 2
    tolerance: Length = 13.0 * u.Angstrom
    w: Length = 3 * u.Angstrom
    mask_width: Length = 20 * u.Angstrom
    cont_width: Length = 70 * u.Angstrom
    fwhm_min: int = 2
    fwhm_max: int = 15


@dataclass
class ConfigOverrides:
    sky: Optional[str] = None
    mask_sky: Optional[str] = None
    line_table: Optional[str] = None
    resolution: Optional[Length] = None
    lines: Union[AllLines, None, List[str]] = None
    fitting: Optional[FittingParametersOverrides] = None


@dataclass
class Config:
    sky: str
    mask_sky: bool
    line_table: str
    resolution: Length
    lines: Union[AllLines, List[str]] = "all"
    fitting: FittingParameters = FittingParameters()
    path: str = "."
    cosmology: Cosmology = Cosmology()

    @property
    def line_list(self):
        table = QTable.read(self.line_table)
        table.sort("wavelength")
        if self.lines == "all":
            return table
        else:
            mask = functools.reduce(
                operator.or_, (table["line"] == line for line in self.lines)
            )
            return table[mask]

    @property
    def sky_list(self):
        if self.sky is None:
            self.mask_sky = False
            return None
        else:
            if self.mask_sky is True:
                return QTable.read(self.sky)
            if self.mask_sky is False:
                return None


@dataclass
class Constants:
    """
    Constants pertaining to cosmology, to fitting procedures and other telescope
    specific parameters. Contains sensible defaults and check for unit type
    correctness.
    """

    path: Optional[str] = None
    cosmology: Optional[Cosmology] = None

    globals: ConfigOverrides = ConfigOverrides()
    setups: Dict[str, ConfigOverrides] = field(default_factory=dict)
    sources: Dict[str, ConfigOverrides] = field(default_factory=dict)

    def __call__(
        self, sample: str, setup_name: str, pointing: str, source_number: int,
    ) -> Config:
        extra = override(self.globals, self.setups.get(setup_name))
        extra = override(
            extra, self.sources.get(f"{sample}.{setup_name}.{pointing}.{source_number}")
        )
        return override(
            Config(**asdict(extra)), {"path": self.path, "cosmology": self.cosmology}
        )


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
