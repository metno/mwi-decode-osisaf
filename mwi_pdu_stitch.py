#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import os, sys
import argparse
import netCDF4
import numpy as np
from shutil import copy2

try:
    c_rows, c_cols = os.popen('stty size', 'r').read().split()
except ValueError as e:
    c_cols = 100

def parse_args():

    parser = argparse.ArgumentParser(description='mwi_dec_stitch',
        formatter_class=lambda prog: argparse.HelpFormatter(
            prog, max_help_position=40, width=int(c_cols)))
    parser.add_argument('-i', '--in', required=False, dest='indir',
                        default=None,
                        help='Input dir if filelist not provided')
    parser.add_argument('-f', '--filelist', required=False, default=None,
                        help='+-separated list of files to merge if input '
                        'directory not set. NOTE +-sign to divide files')
    parser.add_argument('-o', '--out', required=True, dest='outdir',
                        help='Output dir')
    parser.add_argument('-O', '--overwrite', dest='owrite', action='store_true',
                        default=False,
                        help='Overwrite existing outputfile')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        default=False,
                        help='Print more stuff to stdout')
    args = parser.parse_args()

    return args


def mwi_pdu_stitch(flist, outdir, overwrite=False, verbose=False):

    # TODO: Add optional start time and time period to concatenate etc

    # Checking the files
    infiles = []
    tot_scanl = 0
    tot_scanl_h = 0
    for f in flist:
        try:
            outfile
        except NameError as e:
            # Set outfile filename to first valid input filename
            outfile = os.path.join(outdir, os.path.basename(f))
        if ( os.path.isfile(outfile) and not overwrite ):
            print('\n{} already exists, use overwrite option\n'.format(outfile))
            sys.exit(2)
        with netCDF4.Dataset(f, 'r') as nc_in:
            tot_scanl += nc_in.dimensions['n_scanl'].size
            tot_scanl_h += nc_in.dimensions['n_scanl_h'].size
        infiles.append(f)

    #print(infiles)
    lvars = ('lat_l', 'lon_l', 'tb18v', 'tb18h', 'tb23v', 'tb23h', 'tb31v',
             'tb31h', 'scanpos')
    hvars = ('lat_h', 'lon_h', 'tb89v', 'tb89h', 'scanpos_h')

    out_vars = {}
    scanl = 0
    scanl_h = 0
    for f in infiles:
        #print(f)
        with netCDF4.Dataset(f, 'r') as nc_in:
            in_udtime = nc_in.variables['dtime'].units
            try:
                base_time
            except NameError as e:
                base_time = in_udtime
            in_scanl = nc_in.dimensions['n_scanl'].size
            in_scanl_h = nc_in.dimensions['n_scanl_h'].size

            try:
                nc_out
            except NameError as e:
                nc_out = netCDF4.Dataset(outfile, 'w')
                nc_out.createDimension('time', 1)
                nc_out.createDimension('n_scanl', tot_scanl)
                nc_out.createDimension('n_scanp',
                                       nc_in.dimensions['n_scanp'].size)
                nc_out.createDimension('pos', 3)
                nc_out.createDimension('n_scanl_h', tot_scanl_h)
                nc_out.createDimension('n_scanp_h',
                                       nc_in.dimensions['n_scanp_h'].size)

            for lv in lvars:
                try:
                    out_vars[lv][scanl:scanl+in_scanl,:] = nc_in.variables[lv][:]
                except KeyError as e:
                    out_vars[lv] = nc_out.createVariable(
                        lv, nc_in.variables[lv].dtype,
                        nc_in.variables[lv].dimensions,
                        fill_value=nc_in.variables[lv]._FillValue, zlib=True)
                    out_vars[lv].setncatts(
                        {key: nc_in.variables[lv].getncattr(key) \
                        for key in nc_in.variables[lv].ncattrs() \
                        if key != '_FillValue'})
                    out_vars[lv][:in_scanl,:] = nc_in.variables[lv][:]

            for hv in hvars:
                try:
                    out_vars[hv][scanl_h:scanl_h+in_scanl_h,:] = nc_in.variables[hv][:]
                except KeyError as e:
                    out_vars[hv] = nc_out.createVariable(
                        hv, nc_in.variables[hv].dtype,
                        nc_in.variables[hv].dimensions,
                        fill_value=nc_in.variables[hv]._FillValue, zlib=True)
                    out_vars[hv].setncatts(
                        {key: nc_in.variables[hv].getncattr(key) \
                        for key in nc_in.variables[hv].ncattrs() \
                         if key != '_FillValue'})
                    out_vars[hv][:in_scanl_h,:] = nc_in.variables[hv][:]

            try:
                out_vars['time']
            except KeyError as e:
                out_vars['time'] = nc_out.createVariable(
                    'time', nc_in.variables['time'].dtype, ('time'),
                    fill_value=nc_in.variables['time']._FillValue, zlib=True)
                out_vars['time'].setncatts(
                    {key: nc_in.variables['time'].getncattr(key) \
                    for key in nc_in.variables['time'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['time'][:] = nc_in.variables['time'][:]

            try:
                dt = netCDF4.num2date(nc_in.variables['dtime'][:], in_udtime)
                out_vars['dtime'][scanl:scanl+in_scanl] = netCDF4.date2num(
                    dt, base_time)
            except KeyError as e:
                out_vars['dtime'] = nc_out.createVariable(
                    'dtime', nc_in.variables['dtime'].dtype,
                    nc_in.variables['dtime'].dimensions,
                    fill_value=nc_in.variables['dtime']._FillValue, zlib=True)
                out_vars['dtime'].setncatts(
                    {key: nc_in.variables['dtime'].getncattr(key) \
                    for key in nc_in.variables['dtime'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['dtime'][:in_scanl] = nc_in.variables['dtime'][:]

            try:
                dt = netCDF4.num2date(nc_in.variables['dtime_h'][:], in_udtime)
                out_vars['dtime_h'][scanl_h:scanl_h+in_scanl_h] = netCDF4.date2num(dt, base_time)
            except KeyError as e:
                out_vars['dtime_h'] = nc_out.createVariable(
                    'dtime_h', nc_in.variables['dtime_h'].dtype,
                    nc_in.variables['dtime_h'].dimensions,
                    fill_value=nc_in.variables['dtime_h']._FillValue, zlib=True)
                out_vars['dtime_h'].setncatts(
                    {key: nc_in.variables['dtime_h'].getncattr(key) \
                    for key in nc_in.variables['dtime_h'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['dtime_h'][:in_scanl_h] = nc_in.variables['dtime_h'][:]

            try:
                out_vars['n_scanp']
            except KeyError as e:
                out_vars['n_scanp'] = nc_out.createVariable(
                    'n_scanp', nc_in.variables['n_scanp'].dtype,
                    nc_in.variables['n_scanp'].dimensions,
                    fill_value=nc_in.variables['n_scanp']._FillValue, zlib=True)
                out_vars['n_scanp'].setncatts(
                    {key: nc_in.variables['n_scanp'].getncattr(key) \
                    for key in nc_in.variables['n_scanp'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['n_scanp'][:] = nc_in.variables['n_scanp'][:]

            try:
                out_vars['n_scanl'][scanl:scanl+in_scanl] = nc_in.variables['n_scanl'][:] + scanl
            except KeyError as e:
                out_vars['n_scanl'] = nc_out.createVariable(
                    'n_scanl', nc_in.variables['n_scanl'].dtype,
                    nc_in.variables['n_scanl'].dimensions,
                    fill_value=nc_in.variables['n_scanl']._FillValue, zlib=True)
                out_vars['n_scanl'].setncatts(
                    {key: nc_in.variables['n_scanl'].getncattr(key) \
                    for key in nc_in.variables['n_scanl'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['n_scanl'][:in_scanl] = nc_in.variables['n_scanl'][:]

            try:
                out_vars['n_scanp_h']
            except KeyError as e:
                out_vars['n_scanp_h'] = nc_out.createVariable(
                    'n_scanp_h', nc_in.variables['n_scanp_h'].dtype,
                    nc_in.variables['n_scanp_h'].dimensions,
                    fill_value=nc_in.variables['n_scanp_h']._FillValue,
                    zlib=True)
                out_vars['n_scanp_h'].setncatts(
                    {key: nc_in.variables['n_scanp_h'].getncattr(key) \
                    for key in nc_in.variables['n_scanp_h'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['n_scanp_h'][:] = nc_in.variables['n_scanp_h'][:]

            try:
                out_vars['n_scanl_h'][scanl_h:scanl_h+in_scanl_h] = nc_in.variables['n_scanl_h'][:] + scanl_h
            except KeyError as e:
                out_vars['n_scanl_h'] = nc_out.createVariable(
                    'n_scanl_h', nc_in.variables['n_scanl_h'].dtype,
                    nc_in.variables['n_scanl_h'].dimensions,
                    fill_value=nc_in.variables['n_scanl_h']._FillValue,
                    zlib=True)
                out_vars['n_scanl_h'].setncatts(
                    {key: nc_in.variables['n_scanl_h'].getncattr(key) \
                    for key in nc_in.variables['n_scanl_h'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['n_scanl_h'][:in_scanl_h] = nc_in.variables['n_scanl_h'][:]

            try:
                out_vars['scanline'][scanl:scanl+in_scanl,:] = nc_in.variables['scanline'][:] + scanl
            except KeyError as e:
                out_vars['scanline'] = nc_out.createVariable(
                    'scanline', nc_in.variables['scanline'].dtype,
                    nc_in.variables['scanline'].dimensions,
                    fill_value=nc_in.variables['scanline']._FillValue,
                    zlib=True)
                out_vars['scanline'].setncatts(
                    {key: nc_in.variables['scanline'].getncattr(key) \
                    for key in nc_in.variables['scanline'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['scanline'][:in_scanl,:] = nc_in.variables['scanline'][:]

            try:
                out_vars['scanline_h'][scanl_h:scanl_h+in_scanl_h,:] = nc_in.variables['scanline_h'][:] + scanl_h
            except KeyError as e:
                out_vars['scanline_h'] = nc_out.createVariable(
                    'scanline_h', nc_in.variables['scanline_h'].dtype,
                    nc_in.variables['scanline_h'].dimensions,
                    fill_value=nc_in.variables['scanline_h']._FillValue,
                    zlib=True)
                out_vars['scanline_h'].setncatts(
                    {key: nc_in.variables['scanline_h'].getncattr(key) \
                    for key in nc_in.variables['scanline_h'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['scanline_h'][:in_scanl_h,:] = nc_in.variables['scanline_h'][:]

            try:
                out_vars['pos']
            except KeyError as e:
                out_vars['pos'] = nc_out.createVariable(
                    'pos', nc_in.variables['pos'].dtype,
                    nc_in.variables['pos'].dimensions,
                    fill_value=nc_in.variables['pos']._FillValue, zlib=True)
                out_vars['pos'].setncatts(
                    {key: nc_in.variables['pos'].getncattr(key) \
                    for key in nc_in.variables['pos'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['pos'][:] = nc_in.variables['pos'][:]

            try:
                out_vars['sat_xyz'][scanl:scanl+in_scanl,:] = nc_in.variables['sat_xyz'][:]
            except KeyError as e:
                out_vars['sat_xyz'] = nc_out.createVariable(
                    'sat_xyz', nc_in.variables['sat_xyz'].dtype,
                    nc_in.variables['sat_xyz'].dimensions,
                    fill_value=nc_in.variables['sat_xyz']._FillValue, zlib=True)
                out_vars['sat_xyz'].setncatts(
                    {key: nc_in.variables['sat_xyz'].getncattr(key) \
                     for key in nc_in.variables['sat_xyz'].ncattrs() \
                     if key != '_FillValue'})
                out_vars['sat_xyz'][:in_scanl,:] = nc_in.variables['sat_xyz'][:]

            try:
                nc_out.Conventions
                nc_out.end_time_and_date = nc_in.end_time_and_date
                nc_out.from_file = ''
            except AttributeError as e:
                for a in nc_in.ncattrs():
                    nc_out.setncattr(a, nc_in.getncattr(a))

            scanl += in_scanl
            scanl_h += in_scanl_h

    try:
        nc_out.close()
        print('Finished writing', outfile)
    except NameError as e:
        # No files to process?
        pass


if __name__ == '__main__':

    args = parse_args()

    # List of the input files
    if args.filelist is not None:
        if args.indir is not None:
            print("WARNING: Both indir and filelist are set, using the "
                  "filelist and ignoring the indir")

        flist = sorted(args.filelist.split('+'))

    elif args.indir is not None:

        if not os.path.isdir(args.indir):
            print('\nInvalid input dir: {}'.format(args.indir))
            sys.exit(1)

        # Some simple checking
        flist = sorted(os.listdir(args.indir))
        for f in flist:
            if ( not f.startswith('mwi') or f[-3:] != '.nc' ):
                print('Skipping', f)
                continue
        flist = [os.path.join(args.indir, x) for x in flist]

    else:
        print("Either an input directory or a filelist must be provided.")
        sys.exit(1)

    # Checking the output directory
    if not os.path.isdir(args.outdir):
        print('\nInvalid output dir: {}'.format(args.outdir))
        sys.exit(1)

    if os.path.dirname(flist[0]) == os.path.abspath(args.outdir):
        print('\nCannot have same input and output dir\n')
        sys.exit(1)

    mwi_pdu_stitch(flist, args.outdir, overwrite=args.overwrite,
                   verbose=args.verbose)
