#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
from numpy.testing import (
    assert_array_almost_equal,
    assert_array_equal,
)
import scipy.io

import cochlea.zilany2014._zilany2014 as _zilany2014


def test_ihc():

    m = scipy.io.loadmat(
        'data_zilany2014/data_zilany2014.mat',
        squeeze_me=True
    )
    fs = float(m['fs'])
    cf = float(m['cf'])
    sound = m['sound'].astype(float)
    vihc_target = m['vihc']

    vihc = _zilany2014.run_ihc(
        signal=sound,
        cf=cf,
        fs=fs,
        species='cat',
        cohc=1.,
        cihc=1.
    )

    assert_array_almost_equal(
        vihc,
        vihc_target,
        decimal=15
    )




def test_synapse():
    """test_synapse()

    This function has problems, because it's using matlab
    implementation of `resample' and it's difficult to exactly
    reproduce the output.

    During the generation of mat-file the `ffGn' was replaced by
    `zeros' and Matlab's `resample' by self-made downsampling function
    equivalent to python's implementation:

    function resampled = resample_fake(X, P, Q)
        b = fir1(Q, 1/Q);
        a = [1];
        filtered = filtfilt(b, a, X);
        resampled = filtered(1:Q:end);

    """
    m = scipy.io.loadmat(
        'data_zilany2014/data_zilany2014.mat',
        squeeze_me=True
    )
    fs = float(m['fs'])
    cf = float(m['cf'])
    vihc = m['vihc']
    meanrate_target = m['meanrate']

    synout = _zilany2014.run_synapse(
        vihc=vihc,
        fs=fs,
        cf=cf,
        anf_type='hsr',
        powerlaw='approximate',
        ffGn=False
    )
    meanrate = synout / (1 + 0.75e-3*synout)

    assert_array_almost_equal(
        meanrate,
        meanrate_target,
        decimal=10
    )
