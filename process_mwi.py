import os
import sys
import re
import argparse
from argparse import RawDescriptionHelpFormatter
from importlib import import_module, reload
import subprocess
from datetime import datetime
from netCDF4 import Dataset
from satpy import Scene
import numpy as np
import pyresample as pr
import xarray as xr
import pandas as pd
import timeit
import time
from dask.distributed import Client

from mwi_writer import write_osisaf_nc
from decode_mwi import data_select


def parse_args():

    valid_period = ["summer", "winter", "jan", "feb", "mar", "apr", "may",
                    "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    valid_icefmt = ['proc', 'final']

    p = argparse.ArgumentParser('concatenate_paramfiles',
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--infiles', required=True,
                   help='Input MWI files, PLUS-separated (+) because the '
                   'filenames contain commas')
    p.add_argument('-o', '--outname', required=False, default=None,
                   help='Output pathname for the processed file. A default '
                   'name is used if this is not set')
    p.add_argument('-l', '--lmaskpath', required=False, default=None,
                   help='Landmask path if the output should be filtered '
                   'by landmask')

    args = p.parse_args()

    return args


def main():

    tp0 = time.time()

    # Processing the arguments
    args = parse_args()
    # Turn the infiles into a list
    infiles = args.infiles.split('+')

    # Creating an output name if there isn't one, and turning it into a
    # full path if it isn't
    if args.outname is None:
        outname = infiles[0].replace('nc', '_s.nc')
    else:
        outname = args.outname
    if os.path.isfile(outname):
        outname = os.path.join('.', outname)
    outdir = os.path.dirname(outname)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    lmaskpath = args.lmaskpath

    # Call the decoder
    tp1 = time.time()
    mwi_data = data_select(infiles)
    tp2 = time.time()

    # Call the writer
    write_osisaf_nc(mwi_data, outname, lmaskpath=lmaskpath)
    print("Written {}".format(outname))
    tp3 = time.time()

    print("Decoding: {} s".format(tp2 - tp1))
    print("Writing: {} s".format(tp3 - tp2))
    print("Total: {} s".format(tp3 - tp0))

if __name__ == '__main__':

    # Fix from https://stackoverflow.com/questions/72821108/hdf5-warnings-when-accessing-xarray-dataset
    c = Client(n_workers=1, threads_per_worker=1)
    main()
