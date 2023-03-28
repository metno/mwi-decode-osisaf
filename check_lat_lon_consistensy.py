import numpy as np


def check_lat_lon_consistence(lat, lon, ssmis=False, mwri=False, amsr2=False, mwi=False, high_res=False, debug=False):
    limit = 20
    if ssmis:
        limit = 50
    elif mwri:
        limit = 20
    elif amsr2:
        limit = 50
    elif mwi:
        limit = 50
    else:
        print("No instrument given.")
        return

    lat_r = np.radians(lat)
    lon_r = np.radians(lon)

    dlon = lon_r[:, 1:] - lon_r[:, 0:-1]
    dlat = lat_r[:, 1:] - lat_r[:, 0:-1]

    a = np.sin(dlat / 2.)**2 + np.cos(lat_r[:, 0:-1]) * np.cos(lat_r[:, 1:]) * np.sin(dlon / 2.)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1. - a))
    _distance = 6371.0 * c
    distance = np.c_[_distance, _distance[:, -1]]

    ll_q = np.zeros(distance.shape, dtype=np.uint8)
    ll_q[distance > limit] = 1

    if debug:
        if np.any(ll_q):
            print("Some {}latitude longitude are suspicious.".format("high resolution " if high_res else ""))
        else:
            print("All {}latitude longitude scan wise is OK.".format("high resolution " if high_res else ""))

    return ll_q
