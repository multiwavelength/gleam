__author__ = "Andra Stroe"
__version__ = "0.1"

import os, sys
import warnings
import glob
from multiprocessing import Pool

import numpy as np
from astropy.io import fits
import click

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
    (
        extension,
        target,
        inspect,
        fix_center,
        constrain_center,
        bin,
        verbose,
        ignore_sky_lines,
    ) = p
    main.run_main(
        extension,
        target,
        inspect,
        fix_center,
        constrain_center,
        verbose,
        ignore_sky_lines,
        bin,
    )


# Define command line arguments
@click.command()
@click.option("--inspect", is_flag=True)
@click.option("--fix-center", is_flag=True)
@click.option("--constrain-center", is_flag=True)
@click.option("--bin", default=1)
@click.option("--verbose", is_flag=True)
@click.option("--ignore-sky-lines", is_flag=True)
@click.option(
    "--head-path",
    default="/home/andra/Desktop/Keep/Cluster_spectroscopy/line_measurements_mydata_2020_newpipeline_test",
)
@click.option("--max-cpu", default=8, type=int)
def pipeline(
    inspect,
    fix_center,
    constrain_center,
    bin,
    head_path,
    max_cpu,
    verbose,
    ignore_sky_lines,
):
    # Relative paths
    pipeline = f"{head_path}/pipeline"
    data_path = f"{head_path}/allfluxes"
    final_path = f"{head_path}/measurements"
    

    # Define the folder structure such that we run the main line fitting
    # tool in the right place
    unique_sources = []
    for cluster in glob.glob(f"{data_path}/A2254*"):
        for quadrant in glob.glob(f"{cluster}/Q1*"):
            for extension in glob.glob(f"{quadrant}/EXT1*"):
                list_of_targets = (
                    f"{extension}/{os.path.basename(cluster)}_"
                    f"{os.path.basename(quadrant)}_"
                    f"{os.path.basename(extension)}_"
                    f"zinfo.dat"
                )
                targets = rf.read_lof(list_of_targets)
                for target in targets:
                    # if target['Membership']!='member': continue
                    if float(target["SourceNumber"])!= 20:
                        continue
                    # if target["SourceNumber"]!='204': continue
                    # Add all the relevant variable for the running of the
                    # fitting function to a list in the form of a dictionary
                    # This is necessary such that all variable are together
                    # in a single variable for executing the code in a parallel
                    # manner with Pool
                    unique_sources.append(
                        (
                            extension,
                            target,
                            inspect,
                            fix_center,
                            constrain_center,
                            bin,
                            verbose,
                            ignore_sky_lines,
                        )
                    )

    # Set up multithread processing as executing the fitting on different
    # sources is trivially parallelisable
    nproc = 1 if inspect else max_cpu
    with Pool(nproc) as p:
        p.map(run_source, unique_sources)


if __name__ == "__main__":
    pipeline()
