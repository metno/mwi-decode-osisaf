#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def generate_qcvar(tb_dict, channels={'tb19h': 'tb19h', 'tb19v': 'tb19v', 'tb22v': 'tb22v', 'tb37h': 'tb37h',
                                      'tb37v': 'tb37v', 'tb85h': 'tb85h', 'tb85v': 'tb85v'}):
    '''
    Get tb dict with ['19h', '19v', '22v', '37h', '37v', '85h', '85v'] as
    input data

    Return qc_low & qc_high variables
    '''

    # Set tb limits
    tb_lims = {'tb19h': [80.0, 300.0], 'tb19v': [130.0, 310.0],
               'tb22v': [130.0, 310.0],
               'tb37h': [110.0, 300.0], 'tb37v': [130.0, 310.0],
               'tb85h': [110.0, 300.0], 'tb85v': [130.0, 310.0]}

    # Make sure all tb channels are available
    for tb_key in tb_lims.keys():
        if channels[tb_key] not in tb_dict.keys():
            print('tb_dict missing {}'.format(channels[tb_key]))
            return 1  # or something

    # Create QC variables
#    qc_low = np.zeros(tb_dict[channels['tb19v']].shape).astype(np.bool)
#    qc_high = np.zeros(tb_dict[channels['tb85v']].shape).astype(np.bool)
    qc_low = np.zeros(tb_dict[channels['tb19v']].shape).astype(bool)
    qc_high = np.zeros(tb_dict[channels['tb85v']].shape).astype(bool)

    # Set qc_low to True if a channel is out of bounds
    for tb_key in ['tb19h', 'tb19v', 'tb22v', 'tb37h', 'tb37v']:
        qc_low[np.logical_or(tb_dict[channels[tb_key]] < tb_lims[tb_key][0],
                             tb_dict[channels[tb_key]] > tb_lims[tb_key][1])] = True
    # Repeat for qc_high
    for tb_key in ['tb85h', 'tb85v']:
        qc_high[np.logical_or(tb_dict[channels[tb_key]] < tb_lims[tb_key][0],
                              tb_dict[channels[tb_key]] > tb_lims[tb_key][1])] = True

    # Set qc_low to True if a channel polarization difference is out of bounds
    for tb in ['tb19', 'tb37']:
        qc_low[tb_dict[channels[tb + 'v']] - tb_dict[channels[tb + 'h']] < -20.0] = True
    # Repeat for qc_high
    qc_high[tb_dict[channels['tb85v']] - tb_dict[channels['tb85h']] < -20.0] = True

    # Flag entire line in qc_low if more than 10 FOVs are bad
    # qc_low[qc_low.sum(axis=1)>10,:] = True
    # Flag entire line in qc_high if more than 20 FOVs are bad
    # qc_high[qc_high.sum(axis=1)>20,:] = True

    n_scanl, n_scanp = tb_dict[channels['tb19v']].shape
    n_scanl_h, n_scanp_h = tb_dict[channels['tb85v']].shape
    # Flag entire line in qc_low if more than 10% FOVs are bad
    qc_low[qc_low.sum(axis=1) > (n_scanp / 10), :] = True
    # Flag entire line in qc_high if more than 10% FOVs are bad
    qc_high[qc_high.sum(axis=1) > (n_scanp_h / 10), :] = True

    # Remove scanlines that have identical values (check tb19v) for
    # > 10% scanpositions
    for scanl in range(1, n_scanl):
        scanp_limit = (n_scanp - max(tb_dict[channels['tb19v']][scanl].mask.sum(),
                                     tb_dict[channels['tb19v']][scanl - 1].mask.sum())) / 10.0
        if (((tb_dict[channels['tb19v']][scanl, :] -
              tb_dict[channels['tb19v']][scanl - 1, :]) == 0.0).sum() > scanp_limit):
            qc_low[scanl] = True
            qc_low[scanl - 1] = True
    for scanl_h in range(2, n_scanl_h):
        scanp_limit_h = (n_scanp_h - max(tb_dict[channels['tb85v']][scanl].mask.sum(),
                                         tb_dict[channels['tb85v']][scanl - 2].mask.sum())) / 10.0
        if (((tb_dict[channels['tb85v']][scanl_h, :] -
              tb_dict[channels['tb85v']][scanl_h - 2, :]) == 0.0).sum() > scanp_limit_h):
            qc_high[scanl_h] = True
            qc_high[scanl_h - 2] = True

    return (qc_low, qc_high)
