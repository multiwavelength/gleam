r"""
Contains a list of constants and user defined units

"""
__author__ = "Andra Stroe"
__version__ = "0.1"

import yaml
from dataclasses import replace, asdict, is_dataclass, field, dataclass as dat
from typing import Dict, Optional, List, Union, Literal, NamedTuple
import operator
import functools


from pydantic.dataclasses import dataclass
from astropy import units as u
from astropy.table import QTable
from astropy.cosmology import FlatLambdaCDM
from colorama import Fore

AllLines = Literal["all"]
CenterConstraint = Literal["free", "constrained", "fixed"]


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

@dat(frozen=True, eq=True, unsafe_hash=True)
class SourceIdentifier:
    sample: str
    setup: str
    pointing: str
    source: int

    @classmethod
    def __get_validators__(cls):
        yield cls.validate
    
    @classmethod
    def validate(cls, v):
        sample, setup, pointing, source = v.split('.')
        return cls(sample, setup, pointing, int(source))

class Length(Quantity):
    _equivalent_unit = u.m


class Frequency(Quantity):
    _equivalent_unit = (1 / u.s).unit


class Temperature(Quantity):
    _equivalent_unit = u.K


def override(basic, overrides):
    """
    Function to override parameters inside a dictionary or a dataclass. Only
    items set in overrides will replace their corresponding items in basic. 
    Function recursively goes inside each item if they are dict and replaces 
    items inside (element by element).
    Input:
        basic: reference object which we want to override
        overrides: object used to override elements inside basic.
    Return:
        dict with values from basic replaced with values in overrides, where 
        they are set.
    """
    # Is the reference object is empty set it to an empty dict. This is needed
    # to correctly override unset parameters
    if basic is None:
        basic = {}
    # If no overrides are set, just return the reference objct
    if overrides is None:
        return basic
    # Convert inputs to dict.
    if is_dataclass(basic):
        basic = asdict(basic)
    if is_dataclass(overrides):
        overrides = asdict(overrides)
    # Return the override dict, with items replaced with their values from
    # overrides, when they are None.
    return {
        **basic,
        **{
            key: value
            for key, value in overrides.items()
            if value is not None and not isinstance(value, dict)
        },
        **{
            key: override(basic.get(key), value)
            for key, value in overrides.items()
            if isinstance(value, dict)
        },
    }


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
    Class to hold overrides for fitting parameters. The overrides can come from
    setups or individual sources. If no overrides are set, then the fitting
    parameters are set to the defaults in the FittingParameters class.
    """

    SN_limit: Optional[float] = None
    tolerance: Optional[Length] = None
    w: Optional[Length] = None
    mask_width: Optional[Length] = None
    cont_width: Optional[Length] = None
    center: Optional[CenterConstraint] = None


@dataclass
class FittingParameters:
    """
    Paramters needed for the fitting of the lines. Contains sensible default 
    values and checks for correct types.
    """

    SN_limit: float = 2
    tolerance: Length = 26.0 * u.Angstrom
    w: Length = 3 * u.Angstrom
    mask_width: Length = 20 * u.Angstrom
    cont_width: Length = 70 * u.Angstrom
    center: CenterConstraint = "free"


@dataclass
class ConfigOverrides:
    """
    Different configuration overrides. The overrides can come from setups or 
    individual sources. If no overrides are set, then the parameters are set to 
    the defaults and values (if the object is mandatory and no default can be
    set) in the Config class.
    """

    sky: Optional[str] = None
    mask_sky: Optional[str] = None
    line_table: Optional[str] = None
    resolution: Optional[Length] = None
    lines: Union[AllLines, None, List[str]] = None
    fitting: Optional[FittingParametersOverrides] = None


@dataclass
class Config:
    """
    An entire configuration, customizable per telescope and per source.
    """

    line_table: str
    resolution: Length
    sky: Optional[str] = None
    mask_sky: bool = False
    lines: Union[AllLines, List[str]] = "all"
    fitting: FittingParameters = FittingParameters()
    cosmology: Cosmology = Cosmology()

    @property
    def line_list(self):
        """
        Return a QTable containing the lines to be fit for each source.
        Input:
            self.line_table: Fits table on disk
            self.lines: Which lines from the table to select
        Output:
            QTable
        """
        table = QTable.read(self.line_table)
        table.sort("wavelength")
        if self.lines == "all":
            return table
        else:
            lines = {line.strip() for line in self.lines}
            mask = [line.strip() in lines for line in table["line"]]
            return table[mask]

    @property
    def sky_list(self):
        """
        Return a sky table of regions contaminated by sky emission/absorption 
        that should be masked or None if not masking is necessary.
        Input:
            self.line_table: Fits table on disk
            self.mask_sky: Should the wavelengths affected by sky be masked in the
            fitting.
        Output:
            QTable
        """
        # If no sky fits table is set, then no sky lines will be masked
        if self.sky is None:
            if self.mask_sky == True:
                print(
                Fore.YELLOW
                + f"Warning: You want to mask the sky, but no valid sky catalog was set. Sky will not be masked."
            )
            self.mask_sky = False
            return None
        else:
            # If sky fits table exits, decide whether to mask or not based on
            # user preferences.
            if self.mask_sky is True:
                return QTable.read(self.sky)
            if self.mask_sky is False:
                return None


@dataclass
class Constants:
    """
    Constants pertaining to cosmology, to fitting procedures and other telescope
    or source specific parameters. Contains sensible defaults and check for unit 
    type correctness.
    """

    cosmology: Optional[Cosmology] = None

    globals: ConfigOverrides = ConfigOverrides()
    setups: Dict[str, ConfigOverrides] = field(default_factory=dict)
    sources: Dict[SourceIdentifier, ConfigOverrides] = field(default_factory=dict)

    def __call__(
        self, sample: str, setup_name: str, pointing: str, source_number: int,
    ) -> Config:
        # Override the defaults with values set in the used configuration file.
        # First override with any global values, then with telescope specific
        # values and then with any source specific values. Cosmology is fixed
        # for the entire project.
        extra = functools.reduce(
            override,
            [
                {},
                self.globals,
                self.setups.get(setup_name),
                self.sources.get(SourceIdentifier(sample, setup_name, pointing, source_number)),
                {"cosmology": self.cosmology},
            ],
        )

        return Config(**extra)


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
