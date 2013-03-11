#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
from numpy.testing import (
    assert_array_almost_equal,
    assert_array_equal,
)

import cochlea.holmberg2007.traveling_waves as tw



def test_bm_wave():

    data = np.load('data/holmberg2007.npz')

    fs = data['fs']
    sound_scaled = data['sound_scaled']
    xbm_target = data['xbm']
    channel = int(data['channel'])
    cf = tw.real_freq_map[channel]

    xbm = tw.run_bm_wave(
        sound_scaled,
        fs
    )

    assert_array_equal(
        xbm[cf],
        xbm_target
    )





def test_lcr4():

    data = np.load('data/holmberg2007.npz')
    xbm = data['xbm']
    lcr4_target = data['lcr4']
    channel = int(data['channel'])
    cf = tw.real_freq_map[channel]

    lcr4 = tw.run_lcr4(
        xbm=xbm,
        fs=data['fs'],
        cf=cf
    )

    assert_array_equal(
        lcr4,
        lcr4_target
    )



def test_ihcrp():

    data = np.load('data/holmberg2007.npz')
    fs = data['fs']
    lcr4 = data['lcr4']
    ihcrp_target = data['ihcrp']
    channel = int(data['channel'])
    cf = tw.real_freq_map[channel]

    ihcrp = tw.run_ihcrp(
        xbm=lcr4,
        fs=fs,
        cf=cf
    )

    assert_array_equal(
        ihcrp,
        ihcrp_target
    )




def test_ihc_meddis2000():
    """ test_ihc_meddis2000()
    Testing against the DSAM implementation

    """
    data = np.load('data/holmberg2007.npz')

    fs = data['fs']
    ihcrp = data['ihcrp']
    psp_target = data['psp']

    psp = tw.run_ihc_meddis2000(
        ihcrp=ihcrp,
        fs=fs,
        gamma_ca=130,
        beta_ca=400,
        tau_m=1e-4,
        gmax_ca=4.5e-9,
        e_ca=0.066,
        tau_ca=1e-4,
        perm_ca0=0,
        perm_z=2e32,
        power_ca=3,
        replenish_rate_y=10,
        loss_rate_l=2580,
        recovery_rate_r=6580,
        reprocess_rate_x=66.3,
        max_free_pool=8,
    )

    assert_array_almost_equal(
        psp,
        psp_target,
        decimal=17
    )



# def test_an_sg_carney_holmberg2007():

#     data = np.load('data/holmberg2007.npz')

#     fs = data['fs']
#     psp = data['psp']

#     an = tw.run_an_sg_carney_holmberg2007(
#         psp=psp,
#         fs=fs,
#         c0=0.5,
#         c1=0.5,
#         s0=1e-3,
#         s1=12.5e-3,
#         refractory_period=0.75e-3
#     )

#     print(an)
