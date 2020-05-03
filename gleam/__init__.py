__author__ = "Andra Stroe"
__version__ = "0.1"

import os, sys
import warnings
import glob
from multiprocessing import Pool
from functools import reduce
from typing import Tuple

import numpy as np
from astropy.io import fits
from astropy.table import vstack, QTable
import click
from colorama import Fore

import gleam.main
import gleam.read_files as rf

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


class Targets:
    def __init__(self, filter: str) -> None:
        find_masters = glob.glob(filter, recursive=True)
        if not find_masters:
            sys.exit(Fore.RED + "Error! Cannot find any master files.")
        try:
            targets = vstack(
                [rf.read_lof(list_of_targets) for list_of_targets in find_masters]
            )
            targets["key"] = reduce(
                np.char.add,
                (
                    targets["Sample"],
                    ".",
                    targets["Setup"],
                    ".",
                    targets["Pointing"],
                    ".",
                    targets["SourceNumber"].astype(str),
                ),
            )
            targets.add_index("key")
            self._targets: QTable = targets
        except ValueError:
            sys.exit(
                Fore.RED
                + "Error! Cannot stack master files that contain units with those that don't."
            )

    def __getitem__(self, key: Tuple[str, str, str, int]):
        sample, setup, pointing, source = key
        return self._targets.loc[f"{sample}.{setup}.{pointing}.{source}"]


# Define command line arguments
@click.command()
@click.option("--inspect", is_flag=True)
@click.option("--plot", is_flag=True)
@click.option("--fix-center", is_flag=True)
@click.option("--constrain-center", is_flag=True)
@click.option("--bin", default=1)
@click.option("--verbose", is_flag=True)
@click.option("--max-cpu", default=8, type=int)
@click.option("--path", default=".")
@click.option("--spectra")
def pipeline(
    inspect, plot, fix_center, constrain_center, bin, max_cpu, verbose, path, spectra
):
    targets = Targets(f"{path}/**/master*dat")
    if spectra is None:
        spectra = f"{path}/**/spec1d*fits"

    find_spectra = glob.glob(f"{spectra}", recursive=True)
    unique_sources = []

    for spectrum_file in find_spectra:
        _, sample, setup, pointing, source, *_ = os.path.basename(spectrum_file).split(
            "."
        )
        source = int(source)

        try:
            target = targets[sample, setup, pointing, source]
        except KeyError:
            print(
                Fore.RED
                + f"Error! Cannot find source {sample}.{setup}.{pointing}.{source} in any master file."
            )
            continue
        unique_sources.append(
            (
                spectrum_file,
                target,
                inspect,
                plot,
                fix_center,
                constrain_center,
                verbose,
                bin,
            )
        )

    # Set up multithread processing as executing the fitting on different
    # sources is trivially parallelizable
    nproc = 1 if inspect else max_cpu
    with Pool(nproc) as p:
        p.starmap(main.run_main, unique_sources)


if __name__ == "__main__":
    pipeline()
