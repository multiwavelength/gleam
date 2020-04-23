__author__ = "Andra Stroe"
__version__ = "0.1"

import glob

import numpy as np
from astropy.io import fits
from astropy.table import QTable, Column
from astropy import units as u

import constants as c


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

    t.sort("wl_vacuum")
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
    data = np.array(np.genfromtxt(file1, dtype="U8,U20,i,f,f,U10,f,f,U2,U20,U20"))

    # Extract each measurement of interest into a separate Table
    t = QTable()
    t["Setup"] = Column(
        [datum[0] for datum in data],
        dtype="U",
        description="Telescope setup",
    )
    t["Pointing"] = Column(
        [datum[1] for datum in data],
        dtype="U",
        description="Stack type or telescope source",
    )
    t["SourceNumber"] = Column(
        [datum[2] for datum in data], dtype="U", description="SourceNumber"
    )
    t["RA"] = Column(
        [datum[3] for datum in data],
        unit=u.degree,
        dtype="f",
        description="Right Ascension",
    )
    t["DEC"] = Column(
        [datum[4] for datum in data],
        unit=u.degree,
        dtype="f",
        description="Declination",
    )
    t["Cluster"] = Column(
        [datum[5] for datum in data], dtype="U", description="Parent cluster"
    )
    t["Redshift"] = Column(
        [datum[6] for datum in data], dtype="f", description="Redshift"
    )
    t["Confidence"] = Column(
        [datum[7] for datum in data],
        dtype="f",
        description="Confidence as assigned by user. "
        "Scales from 1 (best) to 5 (worst)",
    )
    t["Membership"] = Column([datum[9] for datum in data], dtype="U20")
    t["Type"] = Column([datum[10] for datum in data], dtype="U20", description="Type")
    return t


def naming_convention(data_path, cluster, source_number, type1, mod):
    return "{}/{}.{}_{}.{:03d}.{}".format(
        data_path, type1, cluster, mod, source_number.astype(int), cluster
    )



    