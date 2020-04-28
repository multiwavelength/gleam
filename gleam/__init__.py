__author__ = "Andra Stroe"
__version__ = "0.1"

import os, sys
import warnings
import glob
from multiprocessing import Pool

import numpy as np
from astropy.io import fits
from astropy.table import vstack
import click
from colorama import Fore

import gleam.main
import gleam.read_files as rf
from gleam.constants import a as c

warnings.filterwarnings("ignore")


def run_source(p):
    """
    Convenience function that runs the fitting on a single source. This is 
    needed for the implementation of the multi-threading. The function Pool 
    takes only one parameter, which is contained in the dictionary "p"
    Input:
        p: dictionary that contains the folder name where the data is, the table
           that contains all the details on the source, the line list table that
           will be fit, and a number of click parameters controlling the 
           interactive inspection of the data and the way the center is being
           fit as well as binning the data   
    """
    (extension, target, inspect, fix_center, constrain_center, bin, verbose,) = p
    main.run_main(
        extension, target, inspect, fix_center, constrain_center, verbose, bin,
    )


# Define command line arguments
@click.command()
@click.option("--inspect", is_flag=True)
@click.option("--fix-center", is_flag=True)
@click.option("--constrain-center", is_flag=True)
@click.option("--bin", default=1)
@click.option("--verbose", is_flag=True)
@click.option("--max-cpu", default=8, type=int)
@click.option("--spec-path", default="**/spec1d*fits")
def pipeline(inspect, fix_center, constrain_center, bin, max_cpu, verbose, spec_path):
    find_masters = glob.glob(f"{c.path}/**/master*dat", recursive=True)
    try:
        targets = vstack(
            [rf.read_lof(list_of_targets) for list_of_targets in find_masters]
        )
        for col in [
            "Sample",
            "Setup",
            "Pointing",
            "SourceNumber",
        ]:
            targets.add_index(col)
    except ValueError:
        print(
            Fore.RED
            + "Cannot stack master files that contain units with those that don't."
        )
        sys.exit("Error!")

    find_spectra = glob.glob(f"{c.path}/{spec_path}", recursive=True)
    unique_sources = []
    for spectrum_file in find_spectra:
        _, sample, setup, pointing, source, *_ = os.path.basename(spectrum_file).split(
            "."
        )
        source = int(source)
        target = (
            targets.loc["Sample", sample]
            .loc["Setup", setup]
            .loc["Pointing", pointing]
            .loc["SourceNumber", source]
        )
        unique_sources.append(
            (
                os.path.dirname(spectrum_file),
                target,
                inspect,
                fix_center,
                constrain_center,
                bin,
                verbose,
            )
        )

    # Set up multithread processing as executing the fitting on different
    # sources is trivially parallelizable
    nproc = 1 if inspect else max_cpu
    with Pool(nproc) as p:
        p.map(run_source, unique_sources)


if __name__ == "__main__":
    pipeline()
