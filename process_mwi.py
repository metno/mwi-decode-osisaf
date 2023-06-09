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
from glob import glob

from mwi_writer import write_osisaf_nc
from decode_mwi import data_select


def parse_args():

    valid_period = ["summer", "winter", "jan", "feb", "mar", "apr", "may",
                    "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    valid_icefmt = ['proc', 'final']

    p = argparse.ArgumentParser('concatenate_paramfiles',
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--inputs', required=True,
                   help='Input MWI files, PLUS-separated (+) because the '
                   'filenames contain commas, or else a directory')
    p.add_argument('-p', '--pattern', required=False,
                   default="W_XX-EUMETSAT*nc",
                   help="Input pattern for compiling a list of files")
    p.add_argument('-o', '--outname', required=False, default=None,
                   help='Output pathname for the processed file. A default '
                   'name is used if this is not set. A directory can be '
                   'provided here for the default output name')
    p.add_argument('-l', '--lmaskpath', required=False, default=None,
                   help='Landmask path if the output should be filtered '
                   'by landmask')
    p.add_argument('-c', '--chanlist',
                   default="tb18v,tb18h,tb23v,tb23h,tb31v,tb31h,tb89v,tb89h",
                   help="Select only the required bands in a comma-separated "
                   "list")

    args = p.parse_args()

    return args


def main():

    tp0 = time.time()

    # Processing the arguments
    args = parse_args()
    chanlist = args.chanlist.split(',')

    # Finding a list of the input files:
    if os.path.isdir(args.inputs):
        infiles = sorted(glob(os.path.join(args.inputs, args.pattern)))
    else:
        infiles = args.inputs.split('+')
    outname = args.outname
    lmaskpath = args.lmaskpath

    # Call the decoder
    tp1 = time.time()
    mwi_data = data_select(infiles, chanlist)
    tp2 = time.time()

    # Check if outname is a NetCDF file
    if outname is not None and outname.endswith('.nc'):
        outname = outname
    # If not, create an output name
    else:
        # If outname exists but it not a .nc file, assume it is a directory
        if outname is not None:
            outdir = outname
        else:
            outdir = '.'
        # Set an automatic filename
        datetimes_l = mwi_data['time_start_scan_utc_l']
        datetimes_h = mwi_data['time_start_scan_utc_h']
        swath_beg_datetime_l = (min(dt for dt in datetimes_l if dt.year > 1))
        swath_beg_datetime_h = (min(dt for dt in datetimes_h if dt.year > 1))
        swath_beg_datetime = min(swath_beg_datetime_l,
                                 swath_beg_datetime_h)
        outname = 'mwi_sgb1_n{}_{:%Y%m%d%H%M}_s.nc'.format(
            str(mwi_data['orbit_start']), swath_beg_datetime)
        outname = os.path.join(outdir, outname)
    outdir = os.path.dirname(outname)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Call the writer
    write_osisaf_nc(mwi_data, outname, chanlist, lmaskpath=lmaskpath)
    print("Written {}".format(outname))
    tp3 = time.time()

    print("Decoding: {} s".format(tp2 - tp1))
    print("Writing: {} s".format(tp3 - tp2))
    print("Total: {} s".format(tp3 - tp0))

if __name__ == '__main__':

    # Fix from https://stackoverflow.com/questions/72821108/hdf5-warnings-when-accessing-xarray-dataset
    c = Client(n_workers=1, threads_per_worker=1)
    main()
