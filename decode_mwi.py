import os
import sys
import re
from importlib import import_module, reload
import subprocess
from datetime import datetime
from netCDF4 import Dataset
from satpy import Scene
import numpy as np
import pyresample as pr
import xarray as xr
import pandas as pd
import time

from mwi_info import chan_info, res


def get_band(varname):
    '''Get the band from a variable name that contains e.g. 31v'''
    try:
        bn = re.search('^.*(\d{2}[vh]).*$', varname).group(1)
    except:
        raise ValueError('Variable name {} does not contain a band indicator, '
                         'e.g. 31v'.format(varname))
    for k in chan_info.keys():
        if re.match('tb{}'.format(bn), chan_info[k][2]):
            return(k)
    raise ValueError('Variable name {} does not contain a valid band '
                     'indicator, e.g. 31v'.format(varname))


def llname(band, ll='lat'):
    '''Find the short name and the name in the original data set for the
    latlon groups'''
    for k in chan_info.keys():
        if band == chan_info[k][1]:
            lln = chan_info[k][2].strip('v').strip('h').replace('tb', ll)
            break
    return lln, '{}_pixels_group_{}'.format(ll, band)


def coarsen_along_scanlines(rawdata):
    """ Coarsen the file by averaging along scanlines.
        Scan direction - x
        Flight direction - y
    """

    coarsened_ds = {}
    for item in rawdata.keys():

        # Finding the correct kernel
        if not item.endswith('v') or item.endswith('h'):
            kernel = res[chan_info[get_band('{}v'.format(item))][3]][0]
        else:
            kernel = res[chan_info[get_band(item)][3]][0]
            
            
        # For lat/lon, select values using a stride rather than an averaging.
        # Note: The stride settings are critical to match the length of the
        # coarsened data
        if item[:3] in ('lat', 'lon'):
            lenx = rawdata[item].data.shape[1]
            coarsened_ds[item] = rawdata[item].isel(x=slice(kernel // 2,
                                                    lenx - (kernel // 2),
                                                            kernel))
        # The default is to coarsen by average "kernel" values along scan
        else:
            # If the kernel is an even number, take an average of two windows
            # with one shifted
            if kernel % 2 == 0:
                arrcopy = rawdata[item].copy(deep=True)
                shiftcopy = np.zeros_like(arrcopy.data)
                lena = shiftcopy.shape[1]
                shiftcopy[:, 0: lena - 1] = arrcopy.data[:, 1: lena]
                shiftcopy[:, lena - 1] = arrcopy.data[:, lena - 1]
                arrcopy.data = (arrcopy.data + shiftcopy) / 2.
                arrcopy = arrcopy.coarsen(x=kernel, boundary='trim').mean()
                coarsened_ds[item] = arrcopy
            else:
                coarsened_ds[item] = rawdata[item].coarsen(x=kernel,
                                                        boundary='trim').mean()

    return coarsened_ds


def data_select(infiles):
    """
    Selecting the MWI data to match as closely as possible to SSMIS"
    """

    # Reading in the scene
    reader_kwargs = {'calibration': 'brightness_temperature'}
    scene = Scene(filenames=infiles, reader='mwi_l1b_nc',
                  reader_kwargs=reader_kwargs)

    # Filtering the data
    rawdata = {}
    coarsedata = {}
    seldata = {}
    mwidata = {}

    bandlist = chan_info.keys()
    lllist = list(set([a[1] for a in chan_info.values()]))

    # Read in the uncoarsened data
    # Using only one dask chunk for now
    # https://stackoverflow.com/questions/66935315/disable-xarrays-automatic-use-of-dask-within-a-dask-task
    print("Reading in raw data")
    for llband in lllist:
        for l in ['lat', 'lon']:
            shlname, lname = llname(llband, l)
            scene.load([lname])
            rawdata[shlname] = scene[lname]
            rawdata[shlname] = rawdata[shlname].chunk(chunks={})
    for band in bandlist:
        bandstr = str(band)
        bandname = chan_info[band][2]
        scene.load([bandstr])
        rawdata[bandname] = scene[bandstr]
        rawdata[bandname] = rawdata[bandname].chunk(chunks={})

    # Now coarsening the data
    print("Coarsening the data")
    coarsedata = coarsen_along_scanlines(rawdata)
    
    # Finally collocate the data, using the 31GHz as the base grid (from
    # which both l and h grids are taken)
    print("Collocating the data")
    refband = 5
    refll = chan_info[refband][1]
    kernel_l = res['l'][0]
    # The l grid has already been found for 31GHz
    seldata['lat_l'] = coarsedata[
        chan_info[refband][2].strip('v').strip('h').replace('tb', 'lat')]
    seldata['lon_l'] = coarsedata[
        chan_info[refband][2].strip('v').strip('h').replace('tb', 'lon')]
    # Creating an h grid with half spacing to the l grid. Here the grid
    # starts from "kernel_l // 2" rather than "kernel / 2", so it coincides
    # with the l grid
    kernel_h = res['h'][0]
    lenx = rawdata[chan_info[refband][2].strip('v').strip('h').replace(
        'tb', 'lat')].data.shape[1]
    seldata['lat_h'] = rawdata[
        chan_info[refband][2].strip('v').strip('h').replace('tb', 'lat')
    ].isel(x=slice(kernel_l // 2, lenx - (kernel_h // 2), kernel_h))
    seldata['lon_h'] = rawdata[
        chan_info[refband][2].strip('v').strip('h').replace('tb', 'lon')
    ].isel(x=slice(kernel_l // 2, lenx - (kernel_h // 2), kernel_h))

    # Cycling over the different latlon grids
    for ll in lllist:
        # Find the bands which use this latlon grid
        bandstoproc = []
        for band in bandlist:
            if chan_info[band][1] == ll:
                bandstoproc.append(band)
        # If the latlon grid is the same as the reference grid (i.e. 31V and
        # 31H), then no regridding needed
        if ll == refll:
            for band in bandstoproc:
                bandname = chan_info[band][2]
                print("Do not need to regrid {}".format(bandname))
                seldata[bandname] = coarsedata[bandname]
        # Otherwise find the swaths and regrid
        else:
            for band in bandstoproc:
                bandname = chan_info[band][2]
                print("Regridding {}".format(bandname))
                latname = bandname.strip('v').strip('h').replace('tb', 'lat')
                lonname = bandname.strip('v').strip('h').replace('tb', 'lon')
                swath_def_orig = pr.geometry.SwathDefinition(
                    lons=coarsedata[lonname].data.flatten(),
                    lats=coarsedata[latname].data.flatten())
                if chan_info[band][3] == 'l':
                    targlat = seldata['lat_l'].data.flatten()
                    targlon = seldata['lon_l'].data.flatten()
                    shp = seldata['lat_l'].data.shape
                elif chan_info[band][3] == 'h':
                    targlat = seldata['lat_h'].data.flatten()
                    targlon = seldata['lon_h'].data.flatten()
                    shp = seldata['lat_h'].data.shape
                swath_def_targ = pr.geometry.SwathDefinition(
                    lons=targlon, lats=targlat)
                # Use a radius of influence 100km (tried a smaller ROI of
                # 15km, got many zeros in output)
                # TODO: Check what ROI should be and properly mask zeros
                # TODO - Check epsilon "The distance to a found value is
                # guaranteed to be no further than (1 + eps) times the
                # distance to the correct neighbour. Allowing for uncertanty
                # decreases execution time."
                data_orig_arr = np.array(coarsedata[bandname].data)
                data_orig = data_orig_arr.flatten()
                seldata[bandname] = np.reshape(
                    pr.kd_tree.resample_nearest(swath_def_orig, data_orig,
                                                swath_def_targ,
                                                radius_of_influence=100000,
                                                epsilon=0.5), shp)


    # Reducing along the flight direction. Just in strides for now, but
    # there is a pattern by which it might be useful to use coarsen
    # (weights 1 for the exact point, and 0.5 for the lines forward and back)?
    for item in seldata.keys():
        stride_l = res['l'][1]
        # Finding which stride to use
        if item.endswith('_l'):
            stride = res['l'][1]
        elif item.endswith('_h'):
            stride = res['h'][1]
        elif not item.endswith('v') or item.endswith('h'):
            stride = res[chan_info[get_band('{}v'.format(item))][3]][1]
        else:
            stride = res[chan_info[get_band(item)][3]][1]
        if stride == 1:
            mwidata[item] = seldata[item]
        else:
            leny = seldata[item].data.shape[0]
            # Now this is a numpy array so use a basic slice rather than isel.
            # This doesn't have to mesh with a "coarsen"ed array, so just
            # stride over the whole array.
            mwidata[item] = seldata[item][slice(0, None, stride), :]
            mwidata[item] = seldata[item][slice(0, None, stride), :]

    # Scan start times (1-d array with length flight)
    scene.load(['time_start_scan_utc'])
    leny = scene['time_start_scan_utc'].data.shape[0]
    stride_l = res['l'][1]
    stride_h = res['h'][1]
    tss_utc_l = np.array(scene['time_start_scan_utc'].isel(
        y=slice(0, None, stride_l)))
    mwidata['time_start_scan_utc_l'] = np.array(
        [datetime.utcfromtimestamp(t.astype(datetime) * 1.0e-9)
         for t in tss_utc_l])
    tss_utc_h = np.array(scene['time_start_scan_utc'].isel(
        y=slice(0, None, stride_h)))
    mwidata['time_start_scan_utc_h'] = np.array(
        [datetime.utcfromtimestamp(t.astype(datetime) * 1.0e-9)
         for t in tss_utc_h])

    # Getting a bit more metadata
    # Note - using python netCDF dataset on these files seems OK as it
    # does not open all the data
    orbit_start = []
    processor_version = []
    for infile in infiles:
        dataset = Dataset(infile, 'r')
        orbit_start.append(dataset.__dict__['orbit_start'])
        processor_version.append(dataset['status']['processing'].__dict__['processor_version'])
    mwidata['orbit_start'] = min(orbit_start)
    mwidata['processor_version'] = ','.join(set(processor_version))

    mwidata['from_file'] = ','.join(infiles)

    print("Done with data_select_coarsen")

    return mwidata
