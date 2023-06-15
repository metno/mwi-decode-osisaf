import os
import sys
import numpy as np
import pyresample as pr
from pyorbital import geoloc, astronomy
from matplotlib import pylab as plt

sys.path += [os.path.abspath('../../../../conc_MET/'), ]

# from conc import io_handler

try:
    import numexpr as ne
except ImportError:
    ne = None

# Nominal satellite/instrument altitude from OSCAR
# WARNING: MWI has preliminary value
instr_alt = {'smmr_ni07': 947.0, 'ssmi_f08': 856.0, 'ssmi_f10': 790.0,
             'ssmi_f11': 853.0, 'ssmi_f13': 850.0, 'ssmi_f14': 852.0, 'ssmi_f15': 850.0,
             'ssmi_f16': 850.0, 'ssmi_f17': 850.0, 'ssmi_f18': 850.0, 'ssmi_f19': 850.0,
             'amsr_aq': 705.0, 'amsr_gw1': 700.0, 'mwri_fy3a': 834.0, 'mwri_fy3b': 836.0,
             'mwri_fy3c': 836.0, 'mwri_fy3d': 836.0, 'mwi': 833.0}
# Nominal instrument zenith angle from OSCAR
# WARNING: MWI has filler value
instr_theta = {'smmr_ni07': 50.2, 'ssmi_f08': 53.1, 'ssmi_f10': 53.1,
               'ssmi_f11': 53.1, 'ssmi_f13': 53.1, 'ssmi_f14': 53.1, 'ssmi_f15': 53.1,
               'ssmi_f16': 53.1, 'ssmi_f17': 53.1, 'ssmi_f18': 53.1, 'ssmi_f19': 53.1,
               'amsr_aq': 55.0, 'amsr_gw1': 55.0, 'mwri_fy3a': 53.1, 'mwri_fy3b': 53.1,
               'mwri_fy3c': 53.1, 'mwri_fy3d': 53.1, 'mwi': 55.0}


def get_cartesian_coords(lons, lats, datum='wgs84'):

    orig_shape = None
    if isinstance(lons, np.ndarray) and lons.ndim > 1:
        orig_shape = lons.shape

    lons = lons.ravel()
    lats = lats.ravel()

    coords = np.zeros((lons.size, 3), dtype=lons.dtype)
    deg2rad = lons.dtype.type(np.pi / 180)

    if datum == 'wgs84':
        a = 6378137.
        f = 1. / 298.257223563
        e2 = (2 - f) * f
        if ne:
            N = ne.evaluate("a/sqrt(1.-e2*sin(lats*deg2rad)*sin(lats*deg2rad))")
            Nz = ne.evaluate("N*(1-e2)")
        else:
            N = a / np.sqrt(1. - e2 * np.sin(lats * deg2rad) * np.sin(lats * deg2rad))
            Nz = N * (1 - e2)
    elif datum == 'sphere':
        N = 6370997.0
        Nz = N
    else:
        raise NotImplementedError("Datum {} not implemented".format(datum,))

    if ne:
        coords[:, 0] = ne.evaluate("N*cos(lats*deg2rad)*cos(lons*deg2rad)")
        coords[:, 1] = ne.evaluate("N*cos(lats*deg2rad)*sin(lons*deg2rad)")
        coords[:, 2] = ne.evaluate("Nz*sin(lats*deg2rad)")
    else:
        coords[:, 0] = N * np.cos(lats * deg2rad) * np.cos(lons * deg2rad)
        coords[:, 1] = N * np.cos(lats * deg2rad) * np.sin(lons * deg2rad)
        coords[:, 2] = Nz * np.sin(lats * deg2rad)

    if orig_shape:
        coords = coords.reshape(orig_shape[0], orig_shape[1], 3)

    return coords


def ecef_to_eci(xyz, utc_times):
    """ Convert ECEF (Earth Centered Earth Fixed) positions to ECI (Earth Centered Inertial)
        positions:
            XYZ are cartesian positions in ECEF. Should have shape (...,3)
            UTC_times are UTC_times. Sould have shape (...)

        http://ccar.colorado.edu/ASEN5070/handouts/coordsys.doc

         [X] [C -S 0][X]
         [Y] = [S C 0][Y]
         [Z]eci [0 0 1][Z]ecf

         C and S are cos() and sin() of gmst (Greenwich Meridian Sideral Time)

        Inspired from satellite-js (https://github.com/shashwatak/satellite-js)
    """
    # XYZ and utc_time must have the same shape
    if not xyz.shape[:-1] == utc_times.shape:
        raise ValueError("shape mismatch for XYZ and utc_times (got {} and {})".format(xyz.shape[:-1], utc_times.shape))

    gmst = -1 * astronomy.gmst(utc_times)
    eci = xyz.copy()
    eci[:, 0] = xyz[:, 0] * np.cos(gmst) - xyz[:, 1] * np.sin(gmst)
    eci[:, 1] = xyz[:, 0] * np.sin(gmst) + xyz[:, 1] * np.cos(gmst)
    return eci


def get_satellite_pos(swath_file=None, lons=None, lats=None, swath_nscanpos=None, plot=False, name=None):
    """ Compute the satellite (Cartesian) position
        at each FoV in the swath.
    """

    if swath_file:
        # Read the fields we need
        lons = swath_file.read('lons')
        lats = swath_file.read('lats')
        swath_nscanpos = swath_file.nscanpos
    else:
        swath_nscanpos = lons.shape[1]

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        evl = 10
        maxl = lons.shape[0] - 1
        evs = 5
        print(lats)
        swath_def = pr.geometry.SwathDefinition(lons=lons, lats=lats)
        all_swath_xyz = get_cartesian_coords(lons, lats)
        print(all_swath_xyz.shape)
        if swath_file:
            tb37v = swath_file.read('tb37v')
        else:
            tb37v = np.zeros(lats.shape)
        ax.scatter(all_swath_xyz[:maxl:evl, ::evs, 0], all_swath_xyz[:maxl:evl, ::evs, 1],
                   all_swath_xyz[:maxl:evl, ::evs, 2],
                   '.', c=tb37v[:maxl:evl, ::evs], edgecolor='none')

    # try and see if the files have a sat_xyz variable
    try:
        if swath_file:
            sat_xyz = swath_file.read('sat_xyz')
        else:
            raise KeyError
        sat_xyz_coord = swath_file.readvarattr('sat_xyz', 'coord_system')
        if sat_xyz_coord == 'ecef':
            # the orbit information is given in Earth Centered Earth Fixed (ECEF) coordinate system
            #    so we must transform to Earth Centered Inertial system (and need UTC times for this)
            utc_times = swath_file.read('times')
            if len(utc_times.shape) > 1:
                # handle 2D (scanlines, scanpos) times (we know sat_xyz is 1D (scanlines,1))
                utc_times = utc_times[:, 0]
            sat_xyz = ecef_to_eci(sat_xyz, utc_times)
        elif sat_xyz_coord == 'eci':
            # do nothing
            pass
        else:
            raise ValueError("Got unexpected value {} for :coord_system (variable {})".format(sat_xyz_coord, sat_xyz))

        if plot:
            ax.scatter(sat_xyz[:maxl:evl, 0], sat_xyz[:maxl:evl, 1], sat_xyz[:maxl:evl, 2], '.', c='r')

    except KeyError as k:
        # try to estimate sat_xyz (To Be Checked and Finalized (TL 24/02/2015))
        # modified Amelie Neuville, Thomas Lavergne, 5/06/2018
        # angle changed to 53 ded if (not swath_file)
        # satellite height changed to 850 km  if (not swath_file)
        # centered difference instead of forward difference for the scanpos:
        # scanpos:scanpos+2 becomes scanpos-1:scanpos+2
        # several indexes changes as a consequence of this modifications
        # computation using only the middle of the scan
        if swath_file:
            sat_height = swath_file.altitude * 1000.
            sat_view_angle = 90. - swath_file.theta
        else:
            if name is not None:
                for k in instr_alt.keys():
                    if os.path.basename(name).startswith(k):
                        if k == 'mwi':
                            print("WARNING: MWI altitude and instrument azimuth data is preliminary.")
                        sat_height = instr_alt[k] * 1000.0
                        sat_view_angle = 90.0 - instr_theta[k]
            else:
                sat_height = 850. * 1000.  # 800
                sat_view_angle = 90. - 53.  # 67

        try:
            sat_height
        except NameError:
            print('Satellite type not in altitude/theta dictionaries in geometry.py.')
            raise NameError

        if len(lons.shape) != 2:
            raise ValueError("Expected a 2D swath file (scanlines,scanpos) but longitudes are a 1D variable")
        if lons.shape[1] != swath_nscanpos:
            raise ValueError("The longitudes do not have the expected {} scanpos (got {}). Wrong instrument?".format(
                swath_nscanpos, lons.shape[1]))

        nl = lons.shape[0]
        ns = swath_nscanpos

        # we use only the middle scanpos of the scanline
        scanpositions = np.array([int((ns - 2) / 2)])
        # if several scanpositions used:
        # sat_xyz=np.zeros([len(scanpositions),nl,3])
        iscanpos = 0
        for scanpos in scanpositions:  # in case we want to reverse to using several scan position, keep the loop
            # use pyresample to compute the FoV's Cartesian coordinates
            swath_def = pr.geometry.SwathDefinition(lons=lons[:, scanpos - 1:scanpos + 2],
                                                    lats=lats[:, scanpos - 1:scanpos + 2])
            swath_xyz = swath_def.get_cartesian_coords()  # 3 points along the scan

            # For each scanline, compute the displacement vector from 1 scanpos to the next
            #   this is our approximation to the along-scan direction
            along_scan = swath_xyz[:, 2, :] - swath_xyz[:, 0, :]
            along_scan_length = ((along_scan[:, 0])**2 + (along_scan[:, 1])**2 + (along_scan[:, 2])**2)**0.5
            for d in range(0, 3):
                along_scan[:, d] /= along_scan_length  # unitary 3d vector in direction of the scan

            # the normal to Earth (sphere) is quite obvious to get
            swath_xyz = swath_xyz[:, 1, :]
            sphere_normal = swath_xyz.copy()  # vector normal to the sphere at the chosen scan point
            sphere_normal_length = ((sphere_normal[:, 0])**2 + (sphere_normal[:, 1])**2 + (sphere_normal[:, 2])**2)**0.5
            for d in range(0, 3):
                sphere_normal[:, d] /= sphere_normal_length

            # vector tangent to Earth, pointing towards the platform projectee
            flat_toward = np.cross(sphere_normal, along_scan)  # crossproduct

            # turn the vectors towards the satellite
            vec = np.swapaxes(flat_toward, 0, 1)
            axis = np.swapaxes(along_scan, 0, 1)
            to_sat = geoloc.qrotate(vec, axis, np.deg2rad(-sat_view_angle))
            to_sat = np.swapaxes(to_sat, 0, 1)

            # get satellite position
            O = swath_xyz  # projected position of the observation
            V = to_sat  # direction toward the satellite from this position
            O2 = O[:, 0]**2 + O[:, 1]**2 + O[:, 2]**2
            V2 = V[:, 0]**2 + V[:, 1]**2 + V[:, 2]**2
            OV = O[:, 0] * V[:, 0] + O[:, 1] * V[:, 1] + O[:, 2] * V[:, 2]
            H = sphere_normal_length[0] + sat_height
            Delta = 4 * OV * OV - 4 * V2 * (O2 - H**2)  # Computation of the intersection between a line and a sphere
            L = (-2 * OV + Delta**0.5) / (2. * V2)  # Distance to the sphere
            sat_xyz = O + L[:, None] * V
            # if several scanposition used:
            # sat_xyz[iscanpos,:] = O + L[:,None]*V
            iscanpos = iscanpos + 1
            if plot:
                ax.scatter(sat_xyz[:maxl:evl, 0], sat_xyz[:maxl:evl, 1], sat_xyz[:maxl:evl, 2], '.')

    return sat_xyz


if __name__ == '__main__':

    from pyorbital import version
    print("Use pyorbital {}".format(version.__version__))

    ifile = '/tmp/amsr_aq_200206110308_s.nc'
    # swath_file = io_handler.SwathFactory.get_swath_guess(ifile)
    swath_file = None
    sat_xyz = get_satellite_pos(swath_file, plot=True)
    print(sat_xyz.shape)
