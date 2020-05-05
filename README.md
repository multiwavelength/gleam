# gleam
Galaxy Line Emission &amp; Absorption Modelling

## What is gleam?

**gleam** is a Python package for fitting Gaussian models to emission and absorption lines in large samples of 1D galaxy spectra. **gleam** is tailored to work well without much human interaction on optical and infrared spectra in a wide range of instrument setups and signal-to-noise regimes.

## Features

- Process large numbers of sources in batch mode.
- Jointly fit lines located close together, report upper limits and identify lines without spectral coverage.
- A single configuration file for an entire project, which you can use to customize the fitting for specific telescope/instrument configurations or even single sources.
- Human-readable YAML configuration file, in which units can be specified.
- Plots of an entire spectrum with the line fits, as well as plots for each individual line fitted. 
- Well integrated with Astropy, which enables the use of units and fits tables.
- Uses LMFIT to perform the fitting and report errors on fit parameters.
- Simple installation with pip.

## How to run

**gleam** fits lines in 1D spectra using redshift information from a master file and several other parameters from a central configuration file.

To run **gleam**, the following are needed:
- A set of 1D spectra, in fits table format
- A set of master files, in fits table or ASCII format, where details of each source are listed
- A configuration file, in YAML format, to specify line lists and fitting constraints
- A line catalog, in fits table format
- (Optional) A sky absorption/emission catalog, in fits table format
Details on the input files can be found further down.

To run the **gleam** using the defaults, you can type:
```
gleam
```

**gleam** has a number of optional command line arguments. For details type:
```
gleam --help
```

## Input data and configuration file

### The input spectra
The input spectra should be in fits format, ideally with units in the headers. They should contain 3 columns: the observed wavelength, the flux and the error. In order to identify source across the spectra and the master files, a naming convention needs to be followed:
- Spectrum name: "spec1d.Sample.Setup.Pointing.SourceNumber.fits"

### Master file
The master file contains information about individual sources in the project, such as the setup and pointing they were observed with, the source number to identify them and their redshift. The master file is used to pull information about each source. You can have a single master file or multiple ones, as long as sources are unique between them.

The master file can be in fits format or ASCII format (with commented header), but should contain the following columns:

| Setup | Pointing | SourceNumber | Sample | Redshift |
|-------|----------|--------------|--------|----------| 
| Keck  | P1       | 123          | Cosmos | 1.2303   |

Column descriptions:
- Setup: the telescope, instrument or mode the source was observed with (must match with setup in the spectrum name)
- Pointing: usually multiple pointings or configurations are observed (must match with pointing in the spectrum name),
- SourceNumber: a unique identifier for the source, within a single setup+pointing combination (must match with source number in the spectrum name),
- Sample: parent sample for the source, e.g. if part of a single galaxy cluster or a famous field,
- Redshift: redshift for the source, ideally correct to within 0.0001.

Master file names need to start with the word "master":
- Example master file name: "master.Sample.Setup.Pointing.fits" or "master.Sample.Setup.Pointing.dat"

### Configuration file

The configuration file enables the user to customize the fitting at 4 levels: use the gleam defaults as much as possible, use a global set of parameters for an entire project, specify a number of telescope/instrument specific overrides or even specify overrides for individual sources.

A minimal working configuration example:
```yaml
globals:
  line_table: line_lists/Main_optical_lines.fits
  resolution: 4.4 Angstrom
```

To report luminosities based on the fitted models, **gleam** uses Hubble
constant, Omega0 and the CMB temperature. While these parameters are reasonably
accurate and up to date, your project may require slightly different values. You
can use the `cosmology` section to override one or more of these parameters.
The same cosmology is going to be used consistently across all the spectra
within a project.

Here is another full example of a configuration file demonstrating how you could
define custom cosmological parameters for your project (the values here happen
correspond to the defaults):

```yaml
globals:
  line_table: line_lists/Main_optical_lines.fits
  resolution: 4.4 Angstrom

cosmology:
  H0: 75 km / (Mpc s)
  Om0: 0.3
  Tcmb0: 2.725 K
```

The configuration has a few more fully customizable parameters, related to 
model fitting, line selection and sky absorption masking.

With **gleam**, you can analyze large numbers of spectra in a uniform manner,
even with data taken in different conditions, with different instruments on
different telescopes and for a wide variety of sources. To make it easy to
capture the specifics of each spectrum, **gleam** offers you the possibility to
specify parameters at three different levels.
1. The global level (`globals`) allows you to override the default configuration
   for all the spectra. This is the most coarse level of customization while the
   next levels provide more fine-grained overrides.
2. The setup level offers a way to apply configuration overrides to
   groups of spectra (named `setups`). Each spectrum belongs to exactly one
   setup. What setups mean is entirely up to you. In general, you would use this
   level to capture differences between telescopes or instruments, such as the
   spectral resolution. The configuration parameters specified at this level
   supersede the the global configuration and the built-in defaults.
3. The per-source level allows you to customize the parameters for each and
   every source. While this can be very helpful to account for some particularly
   troublesome cases, it should be used sporadically both due to the associated
   typing burden as well as in the spirit of keeping the results comparable.

The full structure of the configuration file is:

```yaml

globals:
  <global overrides> ...
setups:
  <setup name>:
    <per-setup overrides> ...
  <setup name>:
    <per-setup overrides> ...  
sources:
  <source locator>:
    <per-source overrides> ...
cosmology:
  H0: 75 km / (Mpc s)
  Om0: 0.3
  Tcmb0: 2.725 K
```

The parameters for each spectrum will be computed by stacking the applicable
overrides on top of the default in order: first the global overrides, then the
applicable per-setup overrides (if any) and finally the applicable per-source
overrides.


#### Overrides

Here are the parameters that can be overridden at either the global, setup or
source level. 

#####  Sky bands to be masked. 

There are two parameters that you can use to control whether the fitting should
ignore portions of the spectrum where sky bands may not have been reliably
subtracted. 

First is the path to a catalog which defines the wavelength intervals of sky
bands. The file it points to must be in the fits file format. See below for
exact details of how to create this file.

```yaml
sky: line_lists/Sky_bands.fits
```

The second is a flag that enables or disables the masking of all the sky bands
in the spectrum.

```yaml
mask_sky: True
```

By default (i.e. if no `sky` or `mask_sky` overrides are applied to a source),
there is no sky masking and the entire spectrum is used. In order for masking to
take place, `sky` myst be set appropriately and `mask_sky` must be set to `True`.

Note that the two overrides don't need to be specified at the same level. For
example, you might want to specify `sky` at the `global` level and then just set
`mask_sky` to `True` for individual `sources` (or `setups`) for which the sky
subtraction is inadequate.

##### Emission/absorption lines to fit

By default, it will fit all lines listed in the line table. Otherwise, only
the lines specified under "lines" will be fitted. One needs to use the same 
names for the lines as in the line table.

```yaml
line_table: line_lists/Main_optical_lines.fits # line catalog
lines: # names of line to select from line table
  - Hb 
  - OIII4
  - OIII5  
  - OII
  - Ha
  - NII1
  - SII1
  - SII2
```


##### Fitting parameters 

There are a number of parameters that each affect the way the line models are
fit to the data. To distinguish them visually in the configuration file, they
are grouped under the `fitting` field. However, they can be individually
overridden at any of the three levels (global, per-setup and per-source).

```yaml
fitting:
  # Signal to noise limit for a detection
  SN_limit: 2
  # Tolerance for considering emission lines part of a   group and thus fitting them together
  tolerance: 13.0 Angstrom
  # Range around either side of the Gaussian center to be probed
  w: 3.0 Angstrom
  # Wavelength range used to mask lines
  mask_width: 20.0 Angstrom
  # Range used for selecting continuum left and right of   the source
  cont_width: 70.0 Angstrom
  # Minimum expected FWHM of the line in pixels
  fwhm_min: 2.0
  # Maximum expected FWHM of the line in pixels
  fwhm_max: 15.0
  # Constraints on the center of each gaussian. The options are:
  # - free: (default) the center can be anywhere within the fitting range
  # - constrained: the center must be within a distance of `w` from the expected
  #   position specified in the `line_table`
  # - fixed: the center is fixed to the expected position of the corresponding
  #   line as specified in the `line_table`
  center: constrained
```

### Line catalog

The line catalog contains a list of lines to draw from when fitting.
It should be in the fits format (preferably with units) and contain columns for
the line name, wavelength and (optionally) the LaTeX representation of the line
name (which is only used when plotting).

|line|wavelength|latex|
|-|-|-|
|Ha   | 6564.614  | H$\boldsymbol{\alpha}$ |
|NII1 | 6585.27   | [N{\sc ii}]            |

A subset of the lines can be specified in the configuration file, otherwise
the entire list of lines from the catalog will be used for fitting.

### Sky catalog

In cases where the sky correction is not done perfectly, your data may still be affected by sky absorption or emission. You can specify a list of bands (in observed wavelength units) to avoid. These bands will be masked and disregarded for fitting and treated as if no spectral coverage is available. (i.e. no upper limits will be reported). Masking of the sky can be turned on and off in the config file.

The sky catalog must be in the fits format (preferably with units) and contain
the following columns:

| band | wavelength_min | wavelength_max|
|-|-|-|
| Aband | 7586.0 | 7658.0 |
| Bband | 6864.0 | 6945.0 |

## Installation requirements

- Python 3.8
- pip ^20.0
- (Optional, but recommended) Latex, when plotting.

## How to install
```
pip install git+https://github.com/multiwavelength/gleam
```

