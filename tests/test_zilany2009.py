#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

from unittest import skip

import numpy as np
from numpy.testing import (
    assert_array_almost_equal,
    assert_array_equal,
)
import scipy.io

import cochlea.zilany2009._pycat as _pycat


def test_ihc():

    m = scipy.io.loadmat(
        'data/zilany2009.mat',
        squeeze_me=True
    )
    fs = float(m['Fs'])
    cf = float(m['CF'])
    sound = m['sound'].astype(float)
    vihc_target = m['vihc']

    meout = _pycat.run_middle_ear_filter(
        signal=sound,
        fs=fs
    )
    vihc = _pycat.run_ihc(
        signal=meout,
        cf=cf,
        fs=fs,
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
        'data/zilany2009.mat',
        squeeze_me=True
    )
    fs = float(m['Fs'])
    cf = float(m['CF'])
    vihc = m['vihc']
    synout_target = m['synout']

    synout = _pycat.run_synapse(
        fs=fs,
        vihc=vihc,
        cf=cf,
        anf_type='hsr',
        powerlaw='approximate',
        ffGn=False
    )

    assert_array_almost_equal(
        synout,
        synout_target,
        decimal=9
    )
