__author__ = "Andra Stroe"
__version__ = "0.1"

import glob

import numpy as np
from astropy.io import fits
from astropy.table import QTable, Column
from astropy import units as u

import constants as c


def read_fitsimg(hdulist):
    """ 
    Reads in a fits image with astropy
    Input:
        hdulist: astropy fits image read with astropy fits.open
    Return: 
        Astropy image and header
    """

    try:
        img = hdulist[1].data
        hdr = hdulist[1].header
    except:
        img = hdulist[0].data
        hdr = hdulist[0].header
    return img, hdr


def read_fits_table(data_path):
    """
    Reads in a fits able with astropy. It tries to account for the fact that the 
    data might be in extenstion 1 or 0.
    Input:
        data_path: file location of the fits file containing the information on
                   the lines that will be measured
    Return: 
        Astropy Table
    """

    hdulist = fits.open(data_path)
    try:
        data = hdulist[1].data
    except:
        data = hdulist[0].data
    return data


def read_lol(data_path):
    """
    Reads the list of lines to measure
    Input:
        data_path: file location of the fits file containing the information on
                   the lines that will be measured
    ! Assumes a certain structure for this file! See below 
    Return: 
        Astropy Table of all the information in the line list file: line name, 
        wavelength in vacuum, type (whether it is a single of double Gaussian), 
        (if doublet) separation from the main line and name in Latex format (for 
        convenience, if plotting the lines somewhere)
    """

    data = read_fits_table(data_path)

    t = QTable()
    t["line"] = Column(
        data.field("Line"),
        dtype="U",
        description="Line name, which includes modifier for single or double Gaussian",
    )
    t["wl_vacuum"] = Column(
        data.field("lambda_vacuum"),
        unit=u.Angstrom,
        dtype="f",
        description="Wavelength in vacuum",
    )
    t["type"] = Column(
        data.field("type"), dtype="d", description="Single or double Gaussian"
    )
    t["separation"] = Column(
        data.field("Doublet_sep"),
        unit=u.Angstrom,
        dtype="f",
        description="Separation from the main line, if a doublet",
    )
    t["latex"] = Column(
        data.field("Tex_name"),
        dtype="U",
        description="Line name written in latex format, useful for plotting purposes",
    )

    return t


def read_lof(file1):
    """
    For each cluster, quadrants and extensions, it reads the head file produced 
    by specpro, which contains a list of the sources and their properties
    Input:
        data_path: folder location of head file and all the individual spectra
                   ! assumes the name contains "_zinfo.dat" in order to find it
                   and that there is only 1 such file in each folder
        The format of the head file is the following
        #Column SourceNumber RA DEC Cluster Redshift Confidence Crew Comments
    Return: 
        Astropy Table with measurements of interest: the source number, RA, DEC,
        the parent cluster, redshift (from specpro) and the z confidence
    """
    data = np.array(np.genfromtxt(file1, dtype="U8,i,f,f,U10,f,f,U2,U20,U20"))

    # Extract each measurement of interest into a separate Table
    t = QTable()
    t["Mode"] = Column(
        [datum[0] for datum in data],
        dtype="U",
        description="Stack type or telescope source",
    )
    t["SourceNumber"] = Column(
        [datum[1] for datum in data], dtype="U", description="SourceNumber"
    )
    t["RA"] = Column(
        [datum[2] for datum in data],
        unit=u.degree,
        dtype="f",
        description="Right Ascension",
    )
    t["DEC"] = Column(
        [datum[3] for datum in data],
        unit=u.degree,
        dtype="f",
        description="Declination",
    )
    t["Cluster"] = Column(
        [datum[4] for datum in data], dtype="U", description="Parent cluster"
    )
    t["Redshift"] = Column(
        [datum[5] for datum in data], dtype="f", description="Redshift"
    )
    t["Confidence"] = Column(
        [datum[6] for datum in data],
        dtype="f",
        description="Confidence as assigned by user. "
        "Scales from 1 (best) to 5 (worst)",
    )
    t["Membership"] = Column([datum[8] for datum in data], dtype="U20")
    t["Type"] = Column([datum[9] for datum in data], dtype="U20", description="Type")
    return t


def naming_convention(data_path, cluster, source_number, type1, mod):
    return "{}/{}.{}_{}.{:03d}.{}".format(
        data_path, type1, cluster, mod, source_number.astype(int), cluster
    )


def read_spectrum(data_path, cluster, source_number, mod):
    """
    Read-in the 1D spectrum for 1 source
    Input:
        data_path:  location of the head file and all the individual spectra
                    ! assumes the name has a certain formatting which was set by
                    the naming convention of specpro
                    spec1d.[Cluster]_stack.[SourceNumber].[Cluster].dat
        The format of the 1D spectra file is; units set by hand:
        # LAMBDA FLUX STDEV
        # Angstrom e-16erg/cm2/s/Ang e-16erg/cm2/s/Ang
    Return:
        Astropy Table with wavelength, flux and standard deviation on the flux
        for the source 
    """

    file1 = "{}.dat".format(
        naming_convention(data_path, cluster, source_number, "spec1d", mod)
    )
    spec = np.genfromtxt(file1, dtype="f,f,f", names=True)

    # Extract the measurements into a separate Table
    t = QTable()
    t["wl"] = Column(
        spec["LAMBDA"], unit=u.Angstrom, dtype="f", description="Observed Wavelength"
    )
    t["flux"] = Column(spec["FLUX"], unit=c.fluxunit, dtype="f", description="Flux")
    t["stdev"] = Column(
        spec["STDEV"],
        unit=c.fluxunit,
        dtype="f",
        description="Standard deviation of the 1D flux",
    )

    return t


def read_infofile(data_path, cluster, source_number, mod, target):
    """
    Read info file attached to every source. This info file is in the formal 
    required for specpro. 
    Input:
        data_path:  location of the head file and all the individual spectra
                    ! assumes the name has a certain formatting which was set by
                    the naming convention of specpro
                    info.[Cluster]_stack.[SourceNumber].[Cluster].dat        
    Output:
        RA, DEC
    """
    file1 = "{}.dat".format(
        naming_convention(data_path, cluster, source_number, "info", mod)
    )
    information = np.genfromtxt(file1, dtype="U10, f", unpack=True)
    names = [v[0] for v in information]
    values = [v[1] for v in information]
    for n, v in zip(names, values):
        if n == "RA":
            target["RA"] = v
        if n == "DEC":
            target["DEC"] = v
    return target
