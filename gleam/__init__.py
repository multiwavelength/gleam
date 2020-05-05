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


def find_source_properties(find_spectra, targets):
    for spectrum_file in find_spectra:
        # Get unique sample, setup, pointing and source names
        _, sample, setup, pointing, source, *_ = os.path.basename(spectrum_file).split(
            "."
        )
        source = int(source)

        # Find source is master file
        try:
            target = targets[sample, setup, pointing, source]
        except KeyError:
            print(
                Fore.RED
                + f"Error! Cannot find source {sample}.{setup}.{pointing}.{source} in any master file. Skipping."
            )
            continue
        if isinstance(target, QTable):
            print(type(target))
            print(
                Fore.RED
                + f"Error! Source {sample}.{setup}.{pointing}.{source} appears in multiple master files. Skipping."
            )
            continue
        yield (spectrum_file, target)


# Define command line arguments
@click.command()
@click.option("--path", default=".")
@click.option("--spectra")
@click.option("--config", default="gleamconfig.yaml")
@click.option("--plot", is_flag=True)
@click.option("--inspect", is_flag=True)
@click.option("--verbose", is_flag=True)
@click.option("--bin", default=1)
@click.option("--max-cpu", default=8, type=int)
def pipeline(path, spectra, config, plot, inspect, verbose, bin, max_cpu):
    # Read configuration file
    config = c.read_config(config)

    # Find all the master files as the targets inside them
    targets = Targets(f"{path}/**/master*dat")

    # Find all the spectrum files
    if spectra is None:
        spectra = f"{path}/**/spec1d*fits"
    find_spectra = glob.glob(f"{spectra}", recursive=True)

    # Make a list of all sources with their properties
    unique_sources = (
        (*unique_source, inspect, plot, verbose, bin, config)
        for unique_source in find_source_properties(find_spectra, targets)
    )

    # Set up multithread processing as executing the fitting on different
    # sources is trivially parallelizable
    nproc = 1 if inspect else max_cpu
    with Pool(nproc) as p:
        p.starmap(main.run_main, unique_sources)


if __name__ == "__main__":
    pipeline()
