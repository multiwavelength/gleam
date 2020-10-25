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
import gleam.constants as c

warnings.filterwarnings("ignore")


class Targets:
    """
    Read and stack all metadata files found into a single table.
    """
    def __init__(self, filter: str) -> None:
        find_meta = glob.glob(filter, recursive=True)
        if not find_meta:
            sys.exit(Fore.RED + "Error! Cannot find any metadata files.")
        try:
            targets = vstack(
                [rf.read_lof(list_of_targets) for list_of_targets in find_meta]
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
                + "Error! Cannot stack metadata files that contain units with those that don't."
            )

    def __getitem__(self, key: Tuple[str, str, str, int]):
        """
        Select a unique row of source properties matching the unique identifies
        of each source, ie its sample, setup, pointing and source number. 
        Input:
            tuple with key identifiers of a source
        Return:
            a unique table row with properties of the selected source
        """
        sample, setup, pointing, source = key
        return self._targets.loc[f"{sample}.{setup}.{pointing}.{source}"]


def find_source_properties(find_spectra, targets):
    """
    For a sample of sources, find matches in the metadata file.
    Input:
        find_spectra: globbed list of spectrum file names
        targets: stack of metadata files
    Return:
        unique combination of spectrum and its corresponding properties
    """
    for spectrum_file in find_spectra:
        # Get unique sample, setup, pointing and source names
        _, sample, setup, pointing, source, *_ = os.path.basename(spectrum_file).split(
            "."
        )
        source = int(source)

        # Find source is metadata file
        try:
            target = targets[sample, setup, pointing, source]
        except KeyError:
            print(
                Fore.RED
                + f"Error! Cannot find source {sample}.{setup}.{pointing}.{source} in any metadata file. Skipping."
            )
            continue
        if isinstance(target, QTable):
            print(type(target))
            print(
                Fore.RED
                + f"Error! Source {sample}.{setup}.{pointing}.{source} appears in multiple metadata files. Skipping."
            )
            continue
        yield (spectrum_file, target)


# Define command line arguments
@click.command()
@click.option("--path", default=".", help='Path to recursively look for metadata files and spectra. See --spectra for overrides.')
@click.option("--spectra", help='Filter for spectra file paths. e.g. "./**/spec1d.Cosmos.Keck.P1.*.fits" to select all sources in the Cosmos sample observed with Keck in pointing P1.')
@click.option("--config", default="gleamconfig.yaml", help='Configuration file in YAML format.')
@click.option("--plot", is_flag=True, help='Save plots of spectrum with emission lines fits next to the corresponding spectrum file.')
@click.option("--inspect", is_flag=True, help='Show interactive plots.')
@click.option("--verbose", is_flag=True, help='Print full output from LMFIT.')
@click.option("--bin", default=1, help='Bin the spectrum before fitting.')
@click.option("--nproc", default=8, type=int, help='Number of threads.')
def pipeline(path, spectra, config, plot, inspect, verbose, bin, nproc):
    # Read configuration file
    config = c.read_config(config)

    # Find all the metadata files as the targets inside them
    targets = Targets(f"{path}/**/meta.*")

    # Find all the spectrum files
    if spectra is None:
        spectra = f"{path}/**/spec1d*fits"
    find_spectra = sorted(glob.glob(f"{spectra}", recursive=True))

    # Make a list of all sources with their properties
    unique_sources = (
        (*unique_source, inspect, plot, verbose, bin, config)
        for unique_source in find_source_properties(find_spectra, targets)
    )

    # Set up multithread processing as executing the fitting on different
    # sources is trivially parallelizable
    nproc = 1 if inspect else nproc
    with Pool(nproc) as p:
        p.starmap(gleam.main.run_main, unique_sources)


if __name__ == "__main__":
    pipeline()
