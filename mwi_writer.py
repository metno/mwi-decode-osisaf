import os
import logging
from datetime import date, time, datetime, timedelta
import numpy as np
import numpy.ma as ma

from netCDF4 import Dataset, date2num, num2date

import converter_utils as cu
import ffwmask as ffwm
from geometry import get_satellite_pos
import pyresample as pr
from generate_qcvar import generate_qcvar
try:
    import simplejson as json
except ImportError:
    import json
from check_lat_lon_consistensy import check_lat_lon_consistence

LOG = logging.getLogger(__name__)
FILL_VALUES = {
    'int': 2147483648,
    'long': 2147483648,
    'int32': 2147483648,
    'float': 1.70E+38,
    'double': 1.70E+38,
    'str': 2147483648,
    'f8': -1.0e10,
    'f4': -1.0e10,
    'i2': -32767,
    'i4': 99999999
}
json_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ssmi_sdr_cmsaf_a_b_corr.json')


class BUFR2NetCDFError(Exception):
    pass


def _nc_valid_var_type(var_type):
    return var_type in ['int', 'float', 'str', 'double', 'long']


def netcdf_datatype(type_name):
    if 'float' in type_name:
        return 'd'
    if 'int' in type_name:
        return 'i'
    if 'long' in type_name:
        return 'i'
    if 'string' in type_name:
        return 'b'

    raise BUFR2NetCDFError("Cannot convert %s to NetCDF compatible type" % type_name)


def _create_global_attributes(grp, title, institution, source, history,
                              references, comments):

    LOG.debug("Creating global attributes")

    grp.Conventions = "CF-1.0"
    grp.title = title.encode('latin-1')
    grp.institution = institution.encode('latin-1')
    grp.source = source.encode('latin-1')
    grp.history = history.encode('latin-1')
    grp.references = references.encode('latin-1')
    grp.comments = comments.encode('latin-1')


def timedelta_total_seconds(td):
    """
    Implements timedelta.total_seconds() as available from Py2.7

    :param td: delta time
    :type td: datetime.timedelta object
    :returns: time difference in total number of seconds

    """
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10 ** 6) / 10 ** 6


def write_osisaf_nc(mwidata,
                    outname,
                    lmaskpath=None,
                    plot=False,
                    min_lat=0.,
                    keep_scan=False,
                    no_qc_vals=False):

    outdir = os.path.dirname(outname)
    if not os.path.exists(outdir):
        print("{} Need to create output dir: {}".format(datetime.now(),
                                                        outdir))
        os.makedirs(outdir)

    skip_n90 = False
    dtime_2d = False
    pack_tbs = False
    ssmi_nb_scanpos = 64 # TODO - check

    datetimes_l = mwidata['time_start_scan_utc_l']
    datetimes_h = mwidata['time_start_scan_utc_h']

    osinc_time_unit = "seconds since 1978-01-01 00:00:00"
    osinc_beg_datetime_l = min(dt for dt in datetimes_l if dt.year > 1)
    swath_beg_date_and_time_l = (min(dt for dt in datetimes_l if dt.year > 1))
    swath_end_date_and_time_l = (max(dt for dt in datetimes_l if dt.year > 1))
    osinc_beg_datetime_h = min(dt for dt in datetimes_h if dt.year > 1)
    swath_beg_date_and_time_h = (min(dt for dt in datetimes_h if dt.year > 1))
    swath_end_date_and_time_h = (max(dt for dt in datetimes_h if dt.year > 1))
    osinc_beg_datetime = min(osinc_beg_datetime_l, osinc_beg_datetime_h)
    osinc_beg_time = date2num(osinc_beg_datetime, osinc_time_unit)
    swath_beg_date_and_time = min(swath_beg_date_and_time_l,
                                  swath_beg_date_and_time_h)
    swath_end_date_and_time = max(swath_end_date_and_time_l,
                                  swath_end_date_and_time_h)

    dtime_scn_l = ma.array([timedelta.total_seconds(t - osinc_beg_datetime) for t in datetimes_l if t.year > 1])

    osinc_format = 'reproc'
    dtime_l = dtime_scndtime = dtime_scn_l

    mwidata['lat_l'] = ma.masked_where(mwidata['lat_l'] < -90, mwidata['lat_l'])
    mwidata['lat_l'] = ma.masked_where(mwidata['lat_l'] > 90, mwidata['lat_l'])
    mwidata['lon_l'] = ma.masked_where(mwidata['lon_l'] < -180,
                                       mwidata['lon_l'])
    mwidata['lon_l'] = ma.masked_where(mwidata['lon_l'] > 180, mwidata['lon_l'])
    if not skip_n90:
        mwidata['lat_h'] = ma.masked_where(mwidata['lat_h'] < -90,
                                           mwidata['lat_h'])
        mwidata['lat_h'] = ma.masked_where(mwidata['lat_h'] > 90,
                                           mwidata['lat_h'])
        mwidata['lon_h'] = ma.masked_where(mwidata['lon_h'] < -180,
                                           mwidata['lon_h'])
        mwidata['lon_h'] = ma.masked_where(mwidata['lon_h'] > 180,
                                           mwidata['lon_h'])

    nb_scanline_l = mwidata['lat_l'].shape[0]
    nb_scanpos_l = mwidata['lat_l'].shape[1]
    if not skip_n90:
        nb_scanline_h = mwidata['lat_h'].shape[0]
        nb_scanpos_h = mwidata['lat_h'].shape[1]

    # Build a master QC mask for all data variables: True means 'mask it:
    # not interesting'
    mask_l = np.zeros((nb_scanline_l, nb_scanpos_l)).astype('bool')
    ll_q_l = check_lat_lon_consistence(mwidata['lat_l'], mwidata['lon_l'],
                                       mwi=True)
    if np.any(ll_q_l):
        mask_l[ll_q_l > 0] = True

    if not skip_n90:
        mask_h = np.zeros((nb_scanline_h, nb_scanpos_h)).astype('bool')
        ll_q_h = check_lat_lon_consistence(mwidata['lat_h'], mwidata['lon_h'],
                                           mwi=True)
        if np.any(ll_q_h):
            mask_h[ll_q_h > 0] = True

    # Discard scanlines partly masked.
    # This is done here when the mask is only influenced by lat/lon
    mask_l[mask_l.any(axis=1), :] = True
    mask_h[mask_h.any(axis=1), :] = True

    # TODO - May want some filtering like this after the brightness
    # temperatures are calibrated
#    mwidata['tb18v'] = ma.masked_where(mwidata['tb18v'] < -195,
#                                       mwidata['tb18v'])
#    mwidata['tb18v'] = ma.masked_where(mwidata['tb18v'] > 60,
#                                       mwidata['tb18v'])
#    mwidata['tb18h'] = ma.masked_where(mwidata['tb18h'] < -195,
#                                       mwidata['tb18h'])
#    mwidata['tb18h'] = ma.masked_where(mwidata['tb18h'] > 60,
#                                       mwidata['tb18h'])
#    mwidata['tb31v'] = ma.masked_where(mwidata['tb31v'] < -195,
#                                       mwidata['tb31v'])
#    mwidata['tb31v'] = ma.masked_where(mwidata['tb31v'] > 60,
#                                       mwidata['tb31v'])
#    mwidata['tb31h'] = ma.masked_where(mwidata['tb31h'] < -195,
#                                       mwidata['tb31h'])
#    mwidata['tb31h'] = ma.masked_where(mwidata['tb31h'] > 60,
#                                       mwidata['tb31h'])
#    if not skip_n90:
#        mwidata['tb89v'] = ma.masked_where(mwidata['tb89v'] < -195,
#                                           mwidata['tb89v'])
#        mwidata['tb89v'] = ma.masked_where(mwidata['tb89v'] > 60,
#                                           mwidata['tb89v'])
#        mwidata['tb89h'] = ma.masked_where(mwidata['tb89h'] < -195,
#                                           mwidata['tb89h'])
#        mwidata['tb89h'] = ma.masked_where(mwidata['tb89h'] > 60,
#                                           mwidata['tb89h'])


    celsius_to_kelvin = 273.15

    #tb19v = env_scene.tb19v.astype(float) * 0.01 + celsius_to_kelvin
    #tb19h = env_scene.tb19h.astype(float) * 0.01 + celsius_to_kelvin
    #tb37v = env_scene.tb37v.astype(float) * 0.01 + celsius_to_kelvin
    #tb37h = env_scene.tb37h.astype(float) * 0.01 + celsius_to_kelvin
    #tb22v = env_scene.tb22.astype(float) * 0.01 + celsius_to_kelvin
    #surf_l = np.array(env_scene.surf)

    scanline_l = np.arange(nb_scanline_l)
    scanline_l = np.vstack((scanline_l,) * nb_scanpos_l).T
    scanpos_l = np.arange(nb_scanpos_l)
    scanpos_l = np.vstack((scanpos_l,) * nb_scanline_l)

    # File
    complete = False
    outname_tmp = outname + '.tmp'
    ds = Dataset(outname_tmp, 'w', format='NETCDF4')
    try:
        created = "Created {:%Y-%m-%d}".format(datetime.utcnow())
        _create_global_attributes(ds, 'L1 MWI satellite data',
            'Ocean and Sea Ice Satellite Application Facility (OSI SAF)', '',
            created, '',
            'Brightness Temperatures, discarded fully masked scanlines')

        ds.start_time_and_date = "{:%Y-%m-%dT%H:%M:%SZ}".format(
            swath_beg_date_and_time).encode("latin-1")
        ds.end_time_and_date = "{:%Y-%m-%dT%H:%M:%SZ}".format(
            swath_end_date_and_time).encode("latin-1")
        # TODO - FIX INITIAL FILENAME
        # ds.from_file = os.path.basename(fname)
        # TODO - FIX ORBIT NUMBER
        # ds.orbit_number = filename_data.start_orbit
        ds.satellite = "metopsg-b1"
        ds.scanline_length = nb_scanpos_l
        ds.orbit_start = mwidata['orbit_start']
        ds.processor_version = mwidata['processor_version']
        ds.from_file = mwidata['from_file']

        ds.createDimension('time', 1)
        ds.createDimension('n_scanl', 0)  # Unlimited dimension
        ds.createDimension('n_scanp', nb_scanpos_l)
        ds.createDimension('pos', 3)  # for cartesian positions
        if not skip_n90:
            nb_scanline_h = mwidata['lat_h'].shape[0]
            nb_scanpos_h = mwidata['lat_h'].shape[1]
            ds.createDimension('n_scanl_h', 0)  # Unlimited dimension
            ds.createDimension('n_scanp_h', nb_scanpos_h)
            ds.scanline_length_h = nb_scanpos_h
        else:
            ds.scanline_length_h = "N/A"

        scanline_h = np.arange(nb_scanline_h)
        scanline_h = np.vstack((scanline_h,) * nb_scanpos_h).T
        scanpos_h = np.arange(nb_scanpos_h)
        scanpos_h = np.vstack((scanpos_h,) * nb_scanline_h)

        n_scanl_l = np.arange(nb_scanline_l)
        n_scanp_l = np.arange(nb_scanpos_l)
        n_scanl_h = np.arange(nb_scanline_h)
        n_scanp_h = np.arange(nb_scanpos_h)

        sat_xyz = get_satellite_pos(lons=mwidata['lon_l'],
                                    lats=mwidata['lat_l'],
                                    swath_nscanpos=n_scanp_l, name='mwi')

        dtime_scn_h = ma.array([timedelta.total_seconds(t - osinc_beg_datetime)
                              for t in datetimes_h if t.year > 1])
        dtime_h = dtime_scn_h

        # Create variables with dimensions and attributes.
        # - time
        snc_time = ds.createVariable('time', 'f8', ('time',),
                                     fill_value=cu.getFillValue('f8'))
        snc_time[0] = osinc_beg_time
        snc_time.units = osinc_time_unit
        snc_time.long_name = "reference time of swath data"

        # - scanline
        snc_n_scanl = ds.createVariable('n_scanl', 'i2', ('n_scanl'),
                                        fill_value=cu.getFillValue('i2'))
        snc_n_scanl.units = "1"
        snc_n_scanl.long_name = "Scanline number"

        # - scanpos
        snc_n_scanp = ds.createVariable('n_scanp', 'i2', ('n_scanp'),
                                        fill_value=cu.getFillValue('i2'))
        snc_n_scanp.units = "1"
        snc_n_scanp.long_name = "Position along a scanline"

        # - scanline
        snc_scanline = ds.createVariable('scanline', 'i2', ('n_scanl', 'n_scanp'),
                                         fill_value=cu.getFillValue('i2'))
        snc_scanline.units = "1"
        snc_scanline.long_name = "Scanline number"

        # - scanpos
        snc_scanpos = ds.createVariable('scanpos', 'i2', ('n_scanl', 'n_scanp'),
                                        fill_value=cu.getFillValue('i2'))
        snc_scanpos.units = "1"
        snc_scanpos.long_name = "Position along scanline, low res obs"

        # - pos
        snc_pos = ds.createVariable('pos', 'i2', ('pos',),
                                    fill_value=cu.getFillValue('i2'))
        snc_pos.units = "1"
        snc_pos.long_name = "XYZ dimensions"

        # - lat_l
        snc_lat_l = ds.createVariable('lat_l', 'f4', ('n_scanl', 'n_scanp'),
                                      fill_value=cu.getFillValue('f4'))
        snc_lat_l.units = "degrees_north"
        snc_lat_l.long_name = "latitude low res obs"
        snc_lat_l.datum = "wgs84"

        # - lon_l
        snc_lon_l = ds.createVariable('lon_l', 'f4', ('n_scanl', 'n_scanp'),
                                      fill_value=cu.getFillValue('f4'))
        snc_lon_l.units = "degrees_east"
        snc_lon_l.long_name = "longitude low res obs"
        snc_lon_l.datum = "wgs84"

        # - satellite navigation
        snc_sat_xyz = ds.createVariable('sat_xyz', 'f4', ('n_scanl', 'pos'),
                                        fill_value=cu.getFillValue('f4'))
        snc_sat_xyz.units = "m"
        snc_sat_xyz.long_name = "Cartesian position of spacecraft (at start of scan)"
        snc_sat_xyz.datum = "wgs84"
        snc_sat_xyz.coord_system = 'ecef'  # Earth Centered Earth Fixed

        # - dtime_l
        if dtime_2d:
            snc_dtime = ds.createVariable('dtime', 'i2', ('n_scanl', 'n_scanp'),
                                          fill_value=cu.getFillValue('i2'))
        else:
            snc_dtime = ds.createVariable('dtime', 'i2', ('n_scanl',),
                                          fill_value=cu.getFillValue('i2'))
        try:
            snc_dtime.units = datetime.strftime(num2date(osinc_beg_time,
                                                osinc_time_unit,
                                                only_use_cftime_datetimes=False,
                                                only_use_python_datetimes=True),
                                            'seconds since %Y-%m-%d %H:%M:%S')
        except TypeError:
            snc_dtime.units = datetime.strftime(num2date(osinc_beg_time,
                                                osinc_time_unit),
                                            'seconds since %Y-%m-%d %H:%M:%S')

        snc_dtime.long_name = "time difference from reference time"
        snc_dtime.scale_factor = np.array(1.0).astype('f4')

        if pack_tbs:
            tb_scale, tb_offset = cu.get_scale_offset((90., 300.), (-32700, 32700))
            tb_type = 'i2'
        else:
            tb_type = 'f4'
        tb_fv = cu.getFillValue(tb_type)

        # - tb18v
        snc_tb18v = ds.createVariable('tb18v', tb_type, ('n_scanl', 'n_scanp'),
                                      fill_value=tb_fv)
        snc_tb18v.units = "K"
        snc_tb18v.long_name = "BT 18V GHz"
        snc_tb18v.coordinates = "dtime lat_l lon_l"

        # - tb18h
        snc_tb18h = ds.createVariable('tb18h', tb_type, ('n_scanl', 'n_scanp'),
                                      fill_value=tb_fv)
        snc_tb18h.units = "K"
        snc_tb18h.long_name = "BT 18H GHz"
        snc_tb18h.coordinates = "dtime lat_l lon_l"

        # - tb23v
        snc_tb23v = ds.createVariable('tb23v', tb_type, ('n_scanl', 'n_scanp'),
                                      fill_value=tb_fv)
        snc_tb23v.units = "K"
        snc_tb23v.long_name = "BT 23V GHz"
        snc_tb23v.coordinates = "dtime lat_l lon_l"

        # - tb23h
        snc_tb23h = ds.createVariable('tb23h', tb_type, ('n_scanl', 'n_scanp'),
                                      fill_value=tb_fv)
        snc_tb23h.units = "K"
        snc_tb23h.long_name = "BT 23H GHz"
        snc_tb23h.coordinates = "dtime lat_l lon_l"

        # - tb31v
        snc_tb31v = ds.createVariable('tb31v', tb_type, ('n_scanl', 'n_scanp'),
                                      fill_value=tb_fv)
        snc_tb31v.units = "K"
        snc_tb31v.long_name = "BT 31V GHz"
        snc_tb31v.coordinates = "dtime lat_l lon_l"

        # - tb31h
        snc_tb31h = ds.createVariable('tb31h', tb_type, ('n_scanl', 'n_scanp'),
                                      fill_value=tb_fv)
        snc_tb31h.units = "K"
        snc_tb31h.long_name = "BT 31H GHz"
        snc_tb31h.coordinates = "dtime lat_l lon_l"

        # - tb50v
        #snc_tb50v = ds.createVariable('tb50v', tb_type, ('n_scanl', 'n_scanp'),
        #                              fill_value=tb_fv)
        #snc_tb50v.units = "K"
        #snc_tb50v.long_name = "BT 50V GHz"
        #snc_tb50v.coordinates = "dtime lat_l lon_l"

        # - tb50h
        #snc_tb50h = ds.createVariable('tb50h', tb_type, ('n_scanl', 'n_scanp'),
        #                              fill_value=tb_fv)
        #snc_tb50h.units = "K"
        #snc_tb50h.long_name = "BT 50H GHz"
        #snc_tb50h.coordinates = "dtime lat_l lon_l"

        if not skip_n90:
            # - scanline
            snc_n_scanl_h = ds.createVariable('n_scanl_h', 'i2', ('n_scanl_h'),
                                              fill_value=cu.getFillValue('i2'))
            snc_n_scanl_h.units = "1"
            snc_n_scanl_h.long_name = "Scanline number"

            # - scanpos
            snc_n_scanp_h = ds.createVariable('n_scanp_h', 'i2', ('n_scanp_h'),
                                              fill_value=cu.getFillValue('i2'))
            snc_n_scanp_h.units = "1"
            snc_n_scanp_h.long_name = "Position along a scanline"

            # - scanline
            snc_scanline_h = ds.createVariable('scanline_h', 'i2', ('n_scanl_h', 'n_scanp_h'),
                                               fill_value=cu.getFillValue('i2'))
            snc_scanline_h.units = "1"
            snc_scanline_h.long_name = "Scanline number high res"

            # - scanpos
            snc_scanpos_h = ds.createVariable('scanpos_h', 'i2', ('n_scanl_h', 'n_scanp_h'),
                                              fill_value=cu.getFillValue('i2'))
            snc_scanpos_h.units = "1"
            snc_scanpos_h.long_name = "Position along scanline, high res obs"

            # - lat_h
            snc_lat_h = ds.createVariable('lat_h', 'f4', ('n_scanl_h', 'n_scanp_h'),
                                          fill_value=cu.getFillValue('f4'))
            snc_lat_h.units = "degrees_north"
            snc_lat_h.long_name = "latitude high res obs"
            snc_lat_h.datum = "wgs84"

            # - lon_h
            snc_lon_h = ds.createVariable('lon_h', 'f4', ('n_scanl_h', 'n_scanp_h'),
                                          fill_value=cu.getFillValue('f4'))
            snc_lon_h.units = "degrees_east"
            snc_lon_h.long_name = "longitude high res obs"
            snc_lon_h.datum = "wgs84"

            # - dtime_h
            if dtime_2d:
                snc_dtime_h = ds.createVariable('dtime_h', 'i2', ('n_scanl_h',
                                                                  'n_scanp_h'),
                                            fill_value=cu.getFillValue('i2'))
            else:
                snc_dtime_h = ds.createVariable('dtime_h', 'i2', ('n_scanl_h',),
                                            fill_value=cu.getFillValue('i2'))
            try:
                snc_dtime_h.units = datetime.strftime(num2date(osinc_beg_time,
                                                      osinc_time_unit,
                                                only_use_cftime_datetimes=False,
                                                only_use_python_datetimes=True),
                                            'seconds since %Y-%m-%d %H:%M:%S')
            except TypeError:
                snc_dtime_h.units = datetime.strftime(num2date(osinc_beg_time,
                                                      osinc_time_unit),
                                            'seconds since %Y-%m-%d %H:%M:%S')
            snc_dtime_h.long_name = "time difference from reference time for high res obs"
            snc_dtime_h.scale_factor = np.array(1.0).astype('f4')


            # - tb89v
            snc_tb89v = ds.createVariable('tb89v', tb_type, ('n_scanl_h',
                                                             'n_scanp_h'),
                                          fill_value=tb_fv)
            snc_tb89v.units = "K"
            snc_tb89v.long_name = "BT 89V GHz"
            snc_tb89v.coordinates = "dtime_h lat_h lon_h"

            # - tb89h
            snc_tb89h = ds.createVariable('tb89h', tb_type, ('n_scanl_h',
                                                             'n_scanp_h'),
                                          fill_value=tb_fv)
            snc_tb89h.units = "K"
            snc_tb89h.long_name = "BT 89H GHz"
            snc_tb89h.coordinates = "dtime_h lat_h lon_h"

        # Load the FarFromWaterMask objects if there is a landmask path set
        ffwmasks = []
        if lmaskpath is not None:
            print("Doing ffwm ... ")
            ffwmasks.append(ffwm.MaskFactory.create_mask('ffwm', 'nh',
                                                         datetime.today(),
                                                         lmaskpath=lmaskpath))
            ffwmasks.append(ffwm.MaskFactory.create_mask('ffwm', 'sh',
                                                         datetime.today(),
                                                         lmaskpath=lmaskpath))

        # Mask latitudes below min_lat threshold.
        if min_lat != 0:
            mask_l[abs(mwidata['lat_l']) < min_lat] = True
            mask_h[abs(mwidata['lat_h']) < min_lat] = True

        if len(ffwmasks) != 0:
            # mask on FFWM (Far from water mask)
            swath_def = pr.geometry.SwathDefinition(
                mwidata['lon_l'].reshape(-1),
                mwidata['lat_l'].reshape(-1))
            # print("swath_def: ",swath_def)
            swath_ind = np.arange(mwidata['lon_l'].size, dtype='int')
            # print("swath_ind: ",swath_ind)
            ffwm_keep = np.array([], dtype='int')
            # print("ffwm_keep: ",ffwm_keep)
            try:
                for fmask in ffwmasks:
                    _, ffwm_keep_area = fmask.filter(swath_def, swath_ind)
                    ffwm_keep = np.append(ffwm_keep, ffwm_keep_area)
            except ValueError:
                pass
                # print("Error in ffwm filer: ", ve)

            ffwm_mask = np.array(list(set(swath_ind) - set(ffwm_keep)))
            try:
                ffwm_mask_l = np.unravel_index(ffwm_mask, mask_l.shape)
                mask_l[ffwm_mask_l] = True
            except TypeError:
                # ffwm_mask has shape 0
                mask_l[:, :] = True

            if not skip_n90:
                swath_def_h = pr.geometry.SwathDefinition(
                    mwidata['lon_h'].reshape(-1),
                    mwidata['lat_h'].reshape(-1))
                swath_ind_h = np.arange(mwidata['lon_h'].size, dtype='int')
                ffwm_keep_h = np.array([], dtype='int')
                try:
                    for fmask in ffwmasks:
                        _, ffwm_keep_area_h = fmask.filter(swath_def_h,
                                                           swath_ind_h)
                        ffwm_keep_h = np.append(ffwm_keep_h, ffwm_keep_area_h)
                except ValueError:
                    pass

                ffwm_mask_h = np.array(list(set(swath_ind_h)
                                            - set(ffwm_keep_h)))
                try:
                    ffwm_mask_h = np.unravel_index(ffwm_mask_h, mask_h.shape)
                    mask_h[ffwm_mask_h] = True
                except TypeError:
                    # ffwm_mask has shape 0
                    mask_h[:, :] = True

#        # Calculate quality control for TB data
#        if not no_qc_vals:
#            data = {}
#            data['tb18v'] = tb18v
#            data['tb18h'] = tb18h
#            data['tb37v'] = tb37v
#            data['tb37h'] = tb37h
#            data['tb89v'] = tb89v
#            data['tb89h'] = tb89h
#            qc = generate_qcvar(data)
#            if type(qc) == int:
#                print("Failed to generate qc")
#            else:
#                if qc[0].any():
#                    qc_true = np.where(qc[0] == True)
#                    if mask_lo[qc_true[0], qc_true[1]].all():
#                        print("All flagged qc data already masked")
#                    else:
#                        print("Some value masked in low after qcvar")
#                    mask_lo = np.logical_or(mask_lo, qc[0])
#                if qc[1].any():
#                    qc_high_true = np.where(qc[1] == True)
#                    if mask_hi[qc_high_true[0], qc_high_true[1]].all():
#                        print("All flagged qc high data already masked")
#                    else:
#                        print("Some value masked in high after qcvar")
#                    mask_hi = np.logical_or(mask_hi, qc[1])

        # assign each variable a value
        vals = {}
        ncvar = {}
        #datas = ['lon_l', 'lat_l', 'sat_xyz', 'tb18v', 'tb18h', 'tb23v',
        #         'tb23h', 'tb31v', 'tb31h', 'tb50v', 'tb50h', 'dtime',
        #         'scanline', 'scanpos', 'n_scanl', 'n_scanp']
        datas = ['lon_l', 'lat_l', 'sat_xyz', 'tb18v', 'tb18h', 'tb23v',
                 'tb23h', 'tb31v', 'tb31h', 'dtime',
                 'scanline', 'scanpos', 'n_scanl', 'n_scanp']
        vals['lon_l'] = mwidata['lon_l']
        vals['lat_l'] = mwidata['lat_l']
        vals['sat_xyz'] = sat_xyz
        vals['tb18v'] = mwidata['tb18v']
        vals['tb18h'] = mwidata['tb18h']
        vals['tb23v'] = mwidata['tb23v']
        vals['tb23h'] = mwidata['tb23h']
        vals['tb31v'] = mwidata['tb31v']
        vals['tb31h'] = mwidata['tb31h']
        #vals['tb50v'] = mwidata['tb50v']
        #vals['tb50h'] = mwidata['tb50h']
        vals['dtime'] = dtime_l
        vals['scanline'] = scanline_l
        vals['scanpos'] = scanpos_l
        vals['n_scanl'] = n_scanl_l
        vals['n_scanp'] = n_scanp_l
        ncvar['lon_l'] = snc_lon_l
        ncvar['lat_l'] = snc_lat_l
        ncvar['sat_xyz'] = snc_sat_xyz
        ncvar['tb18v'] = snc_tb18v
        ncvar['tb18h'] = snc_tb18h
        ncvar['tb23v'] = snc_tb23v
        ncvar['tb23h'] = snc_tb23h
        ncvar['tb31v'] = snc_tb31v
        ncvar['tb31h'] = snc_tb31h
        #ncvar['tb50v'] = snc_tb50v
        #ncvar['tb50h'] = snc_tb50h
        ncvar['dtime'] = snc_dtime
        ncvar['scanline'] = snc_scanline
        ncvar['scanpos'] = snc_scanpos
        ncvar['n_scanl'] = snc_n_scanl
        ncvar['n_scanp'] = snc_n_scanp
        if not skip_n90:
            datas += ['lon_h', 'lat_h', 'tb89v', 'tb89h', 'dtime_h',
                      'scanline_h', 'scanpos_h', 'n_scanp_h', 'n_scanl_h']
            vals['lon_h'] = mwidata['lon_h']
            vals['lat_h'] = mwidata['lat_h']
            vals['tb89v'] = mwidata['tb89v']
            vals['tb89h'] = mwidata['tb89h']
            vals['dtime_h'] = dtime_h
            vals['scanline_h'] = scanline_h
            vals['scanpos_h'] = scanpos_h
            vals['n_scanl_h'] = n_scanl_h
            vals['n_scanp_h'] = n_scanp_h
            ncvar['lon_h'] = snc_lon_h
            ncvar['lat_h'] = snc_lat_h
            ncvar['tb89v'] = snc_tb89v
            ncvar['tb89h'] = snc_tb89h
            ncvar['dtime_h'] = snc_dtime_h
            ncvar['scanline_h'] = snc_scanline_h
            ncvar['scanpos_h'] = snc_scanpos_h
            ncvar['n_scanl_h'] = snc_n_scanl_h
            ncvar['n_scanp_h'] = snc_n_scanp_h

        # Discard scanlines fully masked.
        mask_l_subset = mask_l.all(axis=1)
        mask_h_subset = mask_h.all(axis=1)

        all_data_masked = False
        for i, v in enumerate(datas):
            print("Adding to dataset: {}".format(datas[i]))
            if datas[i].startswith('tb') and pack_tbs:
                ncvar[v].set_auto_maskandscale(True)
                ncvar[v].scale_factor = np.float32(tb_scale)
                ncvar[v].add_offset = np.float32(tb_offset)

            if vals[v].shape[0] == nb_scanline_l:
                mask = mask_l
                subset = ~mask_l_subset
            else:
                mask = mask_h
                subset = ~mask_h_subset

            if (v not in ['scanline', 'scanpos', 'sat_xyz', 'lon_l', 'lat_l',
                          'dtime', 'n_scanl', 'n_scanp',
                          'lon_h', 'lat_h', 'scanline_h', 'scanpos_h',
                          'dtime_h', 'n_scanl_h', 'n_scanp_h']):
                try:
                    vals[v].mask = np.logical_or(vals[v].mask, mask)
                except AttributeError:
                    vals[v] = ma.array(vals[v], mask=mask)

                if vals[v].mask.all():
                    print("{} completely masked.".format(datas[i]))
                    all_data_masked = True
                else:
                    all_data_masked = False

            if (keep_scan or v in ['n_scanp', 'n_scanp_h']):
                ncvar[v][:] = vals[v]
            else:
                ncvar[v][:] = vals[v][subset]

        complete = True

    except IOError as ioe:
        print("Some IOError: {}".format(ioe))
        complete = False
    finally:
        ds.close()

    if all_data_masked and complete:
        print("All data masked. Remove this file.")
        os.remove(outname_tmp)
    if len(n_scanl_l) == 0 and complete:
        print("len(n_scanl) == 0. Remove this file.")
        os.remove(outname_tmp)

    if not complete:
        os.remove(outname_temp)

    try:
        os.rename(outname_tmp, outname)
    except OSError:
        pass

    if plot:
        label = {'tb18v': 'Tb 18v (K)', 'tb18h': 'Tb 18h (K)',
                 'tb23v': 'Tb 18v (K)', 'tb23h': 'Tb 18h (K)',
                 'tb31v': 'Tb 31v (K)', 'tb31h': 'Tb 31h (K)',
                 #'tb50v': 'Tb 31v (K)', 'tb50h': 'Tb 31h (K)',
                 'tb89v': 'Tb 89v (K)', 'tb89h': 'Tb 89h (K)'}

        start_time = "{:%Y%m%d%H%M%S}".format(filename_data.start_time)
        print("Start plotting data ... ")
        area_def = pr.utils.load_area('/disk2/pytroll/pyresample/docs/areas.cfg', 'pc_world')
        # Doing l data
        swath_l = pr.geometry.SwathDefinition(lon_l, lat_l)
        #l_data_to_plot = ['tb18v', 'tb18h', 'tb23v', 'tb23h', 'tb31v',
        #                  'tb31h', 'tb50v', 'tb50h']
        l_data_to_plot = ['tb18v', 'tb18h', 'tb23v', 'tb23h', 'tb31v', 'tb31h']
        for l_data in l_data_to_plot:
            print("Doing ", l_data)
            result_l = pr.kd_tree.resample_nearest(swath_l, eval(l_data),
                                                   area_def,
                                                   radius_of_influence=30000,
                                                   fill_value=None)
            pr.plot.save_quicklook(l_data + '-' + start_time + '-quick.png',
                                   area_def, result_l, num_meridians=0,
                                   num_parallels=90, label=label[l_data])

        if not skip_n90:
            swath_h = pr.geometry.SwathDefinition(lon_h, lat_h)
            h_data_to_plot = ['tb89v', 'tb89h']
            for h_data in h_data_to_plot:
                print("Doing ", h_data)
                result_h = pr.kd_tree.resample_nearest(swath_h, eval(h_data),
                                                       area_def,
                                                    radius_of_influence=30000,
                                                    fill_value=None)
                pr.plot.save_quicklook(h_data + '-' + start_time + '-quick.png',
                                       area_def, result_h, num_meridians=0,
                                       num_parallels=90, label=label[h_data])

        print("Done plotting.")
    return True
