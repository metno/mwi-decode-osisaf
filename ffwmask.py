"""

    Simple module for handeling Far From Water masks in the EASE projection. The main entrance
    into this module is the factory method:

    MaskFactory.create_mask(mask_type, hemisphere, timestamp)

"""

import os
import logging

from pyresample import geo_filter, utils, get_area_def

from netCDF4 import Dataset

LOG = logging.getLogger(__name__)


class EASE2Mask(geo_filter.GridFilter):
    """
    Abstract class used as an interface EASE2 mask's with different resolutions.
    The data in the class is retrieved from netcdf files based on the
    hemisphere and timestamp. The location of the masks, the path to the netcdf
    files are defined by each individual subclass.

    :param hemisphere: Either 'nh' for Northern hemispehre or 'sh' for Southern hemisphere
    :type hemisphere: string
    :param timestamp: When the mask is valid
    :type timestamp: datetime

    """

    def __init__(self, hemisphere, timestamp, lmaskpath=None):
        hemisphere = hemisphere.lower()
        pj_str = '+proj=laea +lon_0=0 +ellps=WGS84 +datum=WGS84'
        if hemisphere == 'nh':
            proj_string = pj_str + ' +lat_0=+90'
        elif hemisphere == 'sh':
            proj_string = pj_str + ' +lat_0=-90'
        else:
            raise ValueError('Unknown hemisphere : %s' % hemisphere)

        self.file_path = self.get_file_path(hemisphere, lmaskpath)
        self.mask = self.read_mask(self.file_path, timestamp)

        shape = self.get_shape()
        if self.mask.shape != shape:
            raise ValueError('Unexpected mask shape: ' + str(self.mask.shape))

#        area_def = utils.get_area_def('EASE2 %s' % hemisphere,
        area_def = get_area_def('EASE2 %s' % hemisphere,
                                'EASE2 %s' % hemisphere,
                                'EASE2 %s' % hemisphere,
                                proj_string, shape[0], shape[1],
                                self.get_extent())

        super(EASE2Mask, self).__init__(area_def, self.mask)

    def __str__(self):
        return self.file_path

    def get_file_path(self, hemisphere, lmaskpath):
        raise NotImplementedError('get_file_path method not implemented')

    def read_mask(self, nc_file, timestamp):
        raise NotImplementedError('get_mask method not implemented')

    def get_shape(self):
        raise NotImplementedError('get_shape method not implemented')

    def get_extent(self):
        raise NotImplementedError('get_extent method not implemented')


class EASE2_250Mask(EASE2Mask):
    """
    Abstract base class. Basic EASE2 grid in 25.0 km resolution

    """

    def get_shape(self):
        return (432, 432)

    def get_extent(self):
        return (-5387500., -5387500., +5387500., +5387500.)


class FarFromWaterMask(EASE2_250Mask):
    """
       Defines a mask for discarding geographical positions that are
       far from any water (and are thus not interested for sea ice monitoring)
    """
    def get_file_path(self, hemisphere, lmaskpath):

        fname = 'LandOceanLakeMask_{}_ease2-250.nc'.format(hemisphere)
        if lmaskpath is not None:
            path = os.path.join(lmaskpath, fname)

        else:
            path = os.path.join(os.environ['OSI_ROOT'], 'OSI_HL_Ice/par', fname)
        if not os.path.exists(path):
            raise ValueError("Cannot find mask file %s", path)
        LOG.debug("reading %s from path : %s", self.__class__.__name__, path)
        return path

    def read_mask(self, nc_file, timestamp):
        try:
            r = Dataset(nc_file, 'r')
            mask = r.variables['ffwm'][:]
        except KeyError as k:
            raise ValueError("File %s does not have '%s'" % (nc_file, k,))
        except Exception as e:
            raise ValueError("Cannot open file %s for netCDF reading (%s)" % (nc_file, e,))
        finally:
            r.close()
        return (mask != 0)


class MaskFactory(object):

    @staticmethod
    def create_mask(mask_type, hemisphere, timestamp, lmaskpath=None):
        """
        Factory method. Returns a new mask object based on the input parameters.

        FIXME this method should not be wrappe in a class

        Valid mask types are:

        :ffwm: Mask for discardind positions far from any water

        :param mask_type: The name of the mask
        :type mast_type: string
        :param hemisphere: Either 'nh' for Northern hemispehre or 'sh' for Southern hemisphere
        :type hemisphere: string
        :param timestamp: When the mask is valid
        :type timestamp: datetime

        """
        if mask_type.lower() == 'ffwm':
            mask = FarFromWaterMask(hemisphere, timestamp, lmaskpath=lmaskpath)
            LOG.info("using %s from %s", mask.__class__.__name__, mask)
            return mask

        raise ValueError('Unknown mask type: %s' % mask_type)
