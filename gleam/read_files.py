__author__ = "Andra Stroe"
__version__ = "0.1"

import glob
import sys

import numpy as np
from astropy.io import fits
from astropy.table import QTable, Column
from astropy import units as u
from colorama import Fore

import gleam.constants as c


def read_lof(file1):
    """For each sample, telescope setup and pointing, it reads the master file, 
    which contains a list of the sources and their properties.
    Input:
        file1: master file in ascii or fits format
        The format of the head file is the following
        # Setup Pointing SourceNumber RA DEC Sample Redshift Confidence Initials Membership   Type
    Return: 
        Astropy Table with measurements of interest: the source number, RA, DEC,
        the parent sample, redshift (from specpro) and the z confidence. Throws
        error if the file is not of the right type.
    """
    try:
        table = QTable.read(file1, format="fits")
        return table
    except:
        try:
            table = QTable.read(file1, format="ascii.commented_header")
            return table
        except:
            print(Fore.RED + "Cannot find master redshift file")
            sys.exit("Error!")


def naming_convention(data_path, sample, source_number, setup, pointing, mod):
    """
    Naming convention for files which starts with type of file and is followed
    by details about the source and setup, in order: setup, pointing, source
    number, and parent sample
    """
    return (
        f"{data_path}/{mod}.{setup}.{pointing}.{source_number.astype(int):03d}.{sample}"
    )
