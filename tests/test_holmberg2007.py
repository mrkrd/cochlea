#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
from numpy.testing import assert_array_almost_equal

import cochlea.holmberg2007.traveling_waves as tw



def test_bm_wave():

    data = np.load('data/holmberg2007.npz')

    fs = data['fs']
    sound_scaled = data['sound_scaled']
    xbm_target = data['xbm']
    channel = data['channel']

    xbm = tw.run_bm_wave(
        sound_scaled,
        fs
    )

    assert_array_almost_equal(
        xbm[:,channel],
        xbm_target
    )


def test_ihc_meddis2000():
    """Testing against the DSAM implementation"""

    data = np.load('data/holmberg2007.npz')

    fs = data['fs']
    ihcrp = data['ihcrp']
    psp_target = data['psp']

    psp = tw.run_ihc_meddis2000(
        ihcrp=ihcrp,
        fs=fs,
        gammaCa=130,
        betaCa=400,
        tauCaChan=1e-4,
        GCaMax=4.5e-9,
        CaVrev=0.066,
        tauConcCa=1e-4,
        perm_Ca0=0,
        perm_z=2e32,
        pCa=3,                  # default
        replenishRate_y=10,
        lossRate_l=2580,
        recoveryRate_r=6580,
        reprocessRate_x=66.3,
        maxFreePool_M=8,
        # opmode='spikes'
    )

    assert_array_almost_equal(
        psp,
        psp_target
    )
