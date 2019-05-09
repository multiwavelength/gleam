__author__ = "Andra Stroe"
__version__ = "0.1"

import os, sys
import warnings
import glob
import numpy as np

from astropy.io import fits
import click

import line_fitting as lf
import read_files as rf

warnings.filterwarnings("ignore")

# Define command line arguments
@click.command()
@click.option('--inspect', is_flag=True)
@click.option('--fix-center', is_flag=True)
@click.option('--constrain-center', is_flag=True)
@click.option('--bin', default=1)
@click.option('--head-path', 
       default='/home/andra/Desktop/Keep/Cluster_spectroscopy/ACRes/line_measurements')
@click.option('--lines-table', default='rsvao.fits')#'Main_optical_lines.fits')
def main(inspect, fix_center, constrain_center, bin, head_path, lines_table):
    # Relative paths
    pipeline = f'{head_path}/pipeline'
    data_path = f'{head_path}/allfluxes'
    final_path = f'{head_path}/measurements'
    line_list = rf.read_lol(f'{pipeline}/line_lists/{lines_table}')
    #Define the folder structure such that we run the main line fitting
    #tool in the right place
    clusters = np.genfromtxt(f'{data_path}/clusters_in_my_sample', dtype='U')
    for cluster in glob.glob(f'{data_path}/A115'):#clusters:
        #cluster=f'{data_path}/{cluster}'
        for quadrant in glob.glob(f'{cluster}/Q1*'):
            for extension in glob.glob(f'{quadrant}/EXT1*'):
                list_of_targets = (f"{extension}/{os.path.basename(cluster)}_"
                                  f"{os.path.basename(quadrant)}_"
                                  f"{os.path.basename(extension)}_"
                                  f"zinfo.dat")
                targets = rf.read_lof(list_of_targets)
                for target in targets:#targets[0:1]:
                    if target['Membership']!='member': continue
                    #print(target)
                    #if target["SourceNumber"]!='156': continue
                    print(os.path.basename(cluster), os.path.basename(quadrant), 
                            os.path.basename(extension), 
                            target["SourceNumber"], target["Redshift"], target['Type'])
                    lf.run_main(extension, target, line_list, inspect, 
                                fix_center, constrain_center, bin)
                    # Fix the RA and DEC
                    #target = rf.read_infofile(extension, target["Cluster"], 
                    #                          target["SourceNumber"], 
                    #                          target["Mode"], target)
                    #print('miau')
                # Write out correct list of sources with correct RA and DEC    
                list_of_targets_out = (f"{extension}/{os.path.basename(cluster)}_"
                                  f"{os.path.basename(quadrant)}_"
                                  f"{os.path.basename(extension)}_"
                                  f"zinfo.fits")
                #targets.write(list_of_targets_out, overwrite=True)
                
if __name__ == '__main__':
    main()