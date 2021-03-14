# gleam
Galaxy Line Emission &amp; Absorption Modeling

Citation: [Andra Stroe and Victor-Nicolae Savu 2021 **AJ** 161 158](https://iopscience.iop.org/article/10.3847/1538-3881/abe12a)

[![DOI](https://zenodo.org/badge/230828355.svg)](https://zenodo.org/badge/latestdoi/230828355)

## What is gleam?

**gleam** is a Python package for fitting Gaussian models to emission and absorption lines in large samples of 1D galaxy spectra. **gleam** is tailored to work well without much human interaction on optical and infrared spectra in a wide range of instrument setups and signal-to-noise regimes. **gleam** will create a fits table with Gaussian line measurements, including central wavelength,
width, height and amplitude, as well as estimates for the continuum under the line and the line flux, luminosity, equivalent width and velocity width. **gleam** will also, optionally, make plots of the spectrum with fitted lines overlaid.

## Features

- Process large numbers of sources in batch mode.
- Jointly fit lines located close together, report upper limits and identify lines without spectral coverage.
- A single configuration file for an entire project, which you can use to customize the fitting for specific telescope/instrument configurations or even single sources.
- Human-readable YAML configuration file, in which units can be specified.
- Plots of an entire spectrum with the line fits, as well as plots for each individual line fitted. 
- Well integrated with Astropy, which enables the use of units and fits tables.
- Uses LMFIT to perform the fitting and report errors on fit parameters.
- Simple installation with pip.

## How to run and output

**gleam** fits lines in 1D spectra using redshift information from a metadata file and several other parameters from a central configuration file.

To run **gleam**, the following are needed:
- A set of 1D spectra, in fits table format
- A set of metadata files, in fits table or ASCII format, where details of each source are listed
- A configuration file, in YAML format, to specify line lists and fitting constraints
- A line catalog, in fits table format
- (Optional) A sky absorption/emission catalog, in fits table format

Details on the input files can be found further down.

The outputs of **gleam** include:
- A fits table with line measurements for each source
- (Optional) Plots of the spectrum with overplotted line fits and upper limits.

To run the **gleam** using the defaults, you can type in the terminal:
```
gleam
```

**gleam** has a number of optional command line arguments. For details type:
```
gleam --help
```
 
An example dataset is contained within the git repository. To download it, 
either use the download button or in the terminal:
```
wget https://github.com/multiwavelength/gleam/raw/main/example.tar.gz
```

## Input data and configuration file

### The input spectra
The input spectra should be in fits format, ideally with units in the headers. They should contain 3 columns: the observed wavelength, the flux and the error. In order to identify source across the spectra and the metadata files, a naming convention needs to be followed:
- `spec1d.Sample.Setup.Pointing.SourceNumber.fits`

### Metadata file
The metadata file contains information about individual sources in the project, such as the setup and pointing they were observed with, the source number to identify them and their redshift. The metadata file is used to pull information about each source. You can have a single metadata file or multiple ones, as long as sources are unique between them.

The metadata file can be in fits format or ASCII format (with commented header), but should contain the following columns:

| Setup | Pointing | SourceNumber | Sample | Redshift |
|-------|----------|--------------|--------|----------| 
| Keck  | P1       | 123          | Cosmos | 1.2303   |

Column descriptions:
- **Setup**: the telescope, instrument or mode the source was observed with (must match with setup in the spectrum name)
- **Pointing**: usually multiple pointings or configurations are observed (must match with pointing in the spectrum name),
- **SourceNumber**: a unique identifier for the source, within a single setup+pointing combination (must match with source number in the spectrum name),
- **Sample**: parent sample for the source, e.g. if part of a single galaxy cluster or a famous field (must match with sample in the spectrum name),
- **Redshift**: redshift for the source, ideally correct to within 0.0001.

Metadata files must start with "meta.":
  - `meta.Sample.Setup.Pointing.fits`
  - `meta.Sample.Setup.Pointing.dat`
  - `meta.fits`

### Configuration file

The configuration file enables the user to customize the fitting at 4 levels: use the **gleam** defaults as much as possible, use a global set of parameters for an entire project, specify a number of telescope/instrument specific overrides or even specify overrides for individual sources.

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
   supersede the the global configuration and the built-in defaults. Note that
   the setup name needs to match that in the corresponding sources.
3. The per-source level (named `sources`) allows you to customize the parameters for each and
   every source. While this can be very helpful to account for some particularly
   troublesome cases, it should be used sporadically both due to the associated
   typing burden as well as in the spirit of keeping the results comparable. The
   naming convention of the source should be in line with the input spectrum
   file, without the 'spec1d' and '.fits'. For example: 
   - Example source within sources: "Sample.Setup.Pointing.SourceNumber"


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

#####  Sky bands to be masked 

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
take place, `sky` must be set appropriately and `mask_sky` must be set to `True`.

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
  # Tolerance for considering emission lines part of a group and thus fitting them together
  tolerance: 26.0 Angstrom
  # Range around either side of the Gaussian center to be probed
  w: 3.0 Angstrom
  # Wavelength range used to mask lines
  mask_width: 20.0 Angstrom
  # Range used for selecting continuum left and right of   the source
  cont_width: 70.0 Angstrom
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

## Output

### Line fits tables
For each of the sources in your sample, **gleam** will produce a table with all of the line fits and upper limits (if possible with units derived from the input data). Each line fitted is represented in a separate row, with all the corresponding line fit details contained in different column. The table contains information from the expected wavelength of the line and the redshift of the source, to emission line fit parameters, line fluxes and equivalent widths. 

All of the output files will start with "linefits" and follow the naming convention described above.
- Line fits table: "linefits.Sample.Setup.Pointing.SourceNumber.fits"

An example of the header of such a table and a description of the columns can be found below. Since there are many column, they are listed here in multiple groups.

| line| wavelength |	latex |	z	| zline |	zline_err	| 
|-	|-	|-	|-	|-	|-	|
Ha	| 6564.614 |	H$\boldsymbol{\alpha}$ |	0.26284 |	0.2629235|	3.3654374E-5 |	

| zoffset	| zoffset_err | cont |	cont_err |	wl |	wl_err |	
|-	|-	|-	|-	|-	|-	|
|8.3521445E-5 |	3.3654374E-5 |	0.0012948853	| 1.11991896E-4 |	6565.162	| 0.22092798	|

| height	| height_err	| sigma	| sigma_err	| amplitude |amplitude_err|
|-	|-	|-	|-	|-	|-	|
| 0.017354544	| 0.001047185	| 2.9371586	|0.2215329 |	0.12777047	| 0.008491282	|

| flux |	flux_err |	luminosity |	luminosity_err | EWrest	| EWrest_err|
|-	|-	|-	|-	|-	|-	|
| 0.12777047	| 0.008491282	| 0.27209833	| 0.018082924	| 98.673195	| 10.762495 |

| FWHM|	FWHM_err|	v |	v_err |	detected	|covered |
|-	|-	|-	|-	|-	|-	|
|	6.9164796	| 0.5216701	| 118.18059	| 23.823605 |	true |	true|

A description of each column:
- **line**: Name of the line, from the input line table.
- **wavelength**: Rest wavelength of the line.
- **latex**: Name of the line in Latex format.
- **z**: Redshift of the source, from the input metadata file.
- **zline, zline_err**: Redshift derived from this particular line and its error.
- **zoffset, zoffset_err**: Offset between the systemic source redshift and this line.
- **cont, cont_err**: Continuum estimation around the line and its error.
- **wl, wl_err**: Central wavelength of the line and its error, estimated from the Gaussian fit.
- **height, height_err**: Height of the line and its error, estimated from the Gaussian fit.
- **sigma, sigma_err**: Sigma of the line and its error, estimated from the Gaussian fit.
- **amplitude, amplitude_err**: Amplitude of the line and its error, estimated from the Gaussian fit.
- **flux, flux_err**: Flux of the line and its error.
- **luminosity, luminosity_err**: Luminosity of the line and its error.
- **EWrest, EWrest_err**: Rest-frame equivalent width of the line and its error.
- **FWHM, FWHM_err**: Deconvolved full-width-at-half-maximum width, in wavelength units and its error.
- **v, v_err**: Deconvolved full-width-at-half-maximum velocity width, in wavelength units and its error.
- **detected**: True if line is detected, False if non-detection.
- **covered**: True if line is covered by the spectrum. False is coverage is missing at the location of the lines, i.e. fits or upper limits are not possible.

If the input spectrum has units, the line parameters will also be reported with units.
When the spectral line is not covered by the spectrum, fit values and errors are omitted.
If a line is not detected, **gleam** only reports an upper limit in the amplitude column and omits all other parameters.
The FWHM and the velocity are only reported if the line is spectrally resolved.

### Plots

If plotting is enabled, **gleam** produces two types of figures: a figure showing the entire spectrum with zoom-ins on the emission line fits. The second type of plots are focused on each line fit. Areas masked by sky are shaded gray for clarity.

- Spectrum plot with the fitted lines overlaid: "linefits.Sample.Setup.Pointing.SourceNumber.png"
- Spectrum plot zooming in on Halpha and NII1: "linefits.Sample.Setup.Pointing.SourceNumber.Ha.NII1.png"

*NOTE*: Plotting high quality figures makes **gleam** very slow (a factor of at least 15 slower than without it). Matplotlib with Latex has some memory leak issues, which can cause **gleam** to slowly consume all the memory. I recommend avoiding batch processing
more than 500 sources when also creating plots.

## Installation requirements

- Python 3.8
- pip ^20.0
- (Optional, but recommended) Latex, when plotting.


## How to install
### From pypi
The recommended way to get `gleam` is from the Python Package Index (PyPI). The package is named `astro-gleam` on PyPI:

```
pip install astro-gleam
```

### From source
In rare cases, when you need to install a version that is not published on PyPI, you can install `gleam` directly from the source repository:

```
pip install git+https://github.com/multiwavelength/gleam
```

You can learn about what options are available when installing from source by reading the [official documentation](https://packaging.python.org/tutorials/installing-packages/#installing-from-vcs).

## How to cite **gleam**

If you use **gleam** in your published projects or as a dependency in your code, please include a citation to the companion paper in the Astronomical Journal, as well as a citation to the repository through Zenodo:

Citation: [Andra Stroe and Victor-Nicolae Savu 2021 **AJ** 161 158](https://iopscience.iop.org/article/10.3847/1538-3881/abe12a)

[![DOI](https://zenodo.org/badge/230828355.svg)](https://zenodo.org/badge/latestdoi/230828355)
