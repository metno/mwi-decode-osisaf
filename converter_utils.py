#!/usr/bin/env python

"""

.. module:: converter_utils
   :platform: Unix
   :synopsis: Collection of routines used by the PMW converter fuctions

.. moduleauthor:: Thomas Lavergne <t.lavergne@met.no>, Kristian Rune Larsen <krl@dmi.dk>


"""

import os
import logging

import numpy as np

LOG = logging.getLogger(__name__)


def touch(fname, times=None):
    with file(fname, 'a'):
        os.utime(fname, times)


def timedelta_total_seconds(td):
    """
    Implements timedelta.total_seconds() as available from Py2.7

    :param td: delta time
    :type td: datetime.timedelta object
    :returns: time difference in total number of seconds

    """
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6


undefNC = {'f8': -1.0e10,
           'f4': -1.0e10,
           'i2': -32767,
           'i4': 99999999}


def getFillValue(v):
    """
    Returns an approproate fill value for this datatype. The fill value must be
    set directly in netcdf variable creation using the fill_value keywork

    :param v: string representation of type
    :type v: string

    """
    return undefNC[v]


def validate_var(name, vals):
    """
    Validate data range for a CMSAF SSM/I variable

    The min/max for the TBs are taken from the CMSAF ATBD v1.3 page 28

    :param name: name of the variable
    :type name: string
    :param vals: array of values
    :type vals: numpy.ndarray

    """
    if name[:3] == 'lat':
        return (vals.min() >= -90. and vals.max() <= +90.)
    elif name[:3] == 'lon':
        return (vals.min() >= -180. and vals.max() <= +180.)
    elif name == 'tb19v' or name == 'tb22v' or name == 'tb37v' or name == 'tb90v':
        return (vals.min() >= 130. and vals.max() <= 310.)
    elif name == 'tb19h':
        return (vals.min() >= 80. and vals.max() <= 300.)
    elif name == 'tb37h' or name == 'tb90h':
        return (vals.min() >= 110. and vals.max() <= 300.)
    else:
        raise KeyError("Do not know how to validate %s" % name)


def get_scale_offset(unpacked_range, packed_range):
    """
       Compute netCDF scale_factor and add_offset for packing
       a range of values (typically float) to a range of
       packed values (typically short).

       unpacked_value = packed_value * scale_factor + add_offset

       note: the returned scale_factor and add_offset are float64 and
          should be converted to appropriate types
    """
    # validate input
    if len(unpacked_range) != 2 or len(packed_range) != 2:
        raise ValueError("Expects two ranges of values: unpacked:(a,b) and packed:(c,d)")
    for var in (unpacked_range, packed_range,):
        if (var[-1] <= var[0]):
            raise ValueError("Expects ranges to ordered (min,max)")
    # compute and return
    scale_factor = np.float64(unpacked_range[1] - unpacked_range[0]) / np.float64(packed_range[1] - packed_range[0])
    add_offset = np.float64(unpacked_range[1]) - packed_range[1] * scale_factor
    return scale_factor, add_offset


def DayRange(first_day, last_day):
    """
       Generator for iterating every day between two dates (both included)
    """
    from datetime import timedelta
    for n in range(int((last_day - first_day).days) + 1):
        yield first_day + timedelta(days=n)


def CleanDir(Path, ignore=['.svn', '.CVS']):
    """
       Wipe-out/Clean everything under a top-level path (which is preserved).
    """
    import glob
    import shutil
    if not os.path.isdir(Path):
        raise ValueError("Path given to CleanDir() is not a directory (%s)" % (Path,))

    # wipe-out test_data/output area
    rubbish = glob.glob(Path + '/*')
    for elem in rubbish:
        if os.path.basename(elem) in ignore:
            next
        if os.path.isdir(elem):
            shutil.rmtree(elem)
        else:
            os.remove(elem)


def AlternativeName(fname):
    """
        Find an existing alternative (zipped) name for input file
           if it exists
    """
    if os.path.exists(fname.replace('F0801GL.nc', 'F0802GL.nc')):
        return fname.replace('F0801GL.nc', 'F0802GL.nc')
    if os.path.exists(fname):
        return fname
    for ext in ['.bz2', '.gz']:
        zipped_fname = fname + ext
        if os.path.exists(zipped_fname):
            return zipped_fname
    return 'none'


def Unzip(fname, DestDir=None):
    """
        Apply appropriate unzip, if necessary
        Unzipped file is created in a specified destination directory,
           or a default temp-name if no info is given
    """
    from tempfile import mkdtemp
    for ext, prog in zip(['.bz2', '.gz'], ['bunzip2 -c', 'gunzip -c']):
        if fname.endswith(ext):
            unzipped_fname = fname.rstrip(ext)
            if os.path.exists(unzipped_fname):
                return unzipped_fname, 'none'
            if DestDir is None:
                DestDir = mkdtemp()
            unzipped_fname = os.path.join(DestDir, os.path.basename(unzipped_fname))
            print("Unzip (%s) file %s to %s" % (prog, fname, DestDir))
            status = os.system("%s %s > %s" % (prog, fname, unzipped_fname))
            if status != 0:
                raise ValueError("unable to use %s on %s" % (prog, fname,))
            return unzipped_fname, ext
    return fname, 'none'


def Zip(fname, ext='.bz2'):
    """
        Apply appropriate zipping
    """
    exts = ['.bz2', '.gz']
    zips = ['bzip2', 'gzip']
    try:
        prog = zips[exts.index(ext)]
    except ValueError:
        raise ValueError("Do not know how to zip to %s" % ext)

    zipped_fname = fname + ext
    if os.path.exists(zipped_fname):
        print("Zipped version (%s) of file %s exists already. Do nothing." % (ext, fname,))
        return zipped_fname
    print("Zip file %s to %s" % (fname, ext,))
    status = os.system("%s %s" % (prog, fname,))
    if status != 0:
        raise ValueError("unable to use %s on %s" % (prog, fname,))
    return zipped_fname


def to_bin_str(v):
    """
    @brief return the binary representation of v as a string
    @param v an integer value
    @return the binary representation of v as a string
    """
    if v:
        return to_bin_str(v >> 1) + str(v & 1)
    else:
        return ""


def from_bin_str(s):
    bmsk = 0
    for i, m in enumerate(s[::-1]):
        bmsk += int(m) * 2**i
    return bmsk
