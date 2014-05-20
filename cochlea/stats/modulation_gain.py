#!/usr/bin/env python

"""Modulation gain of inner ear models.

"""

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd

import mrlib as mr
import mrlib.thorns as th
import mrlib.waves as wv

from . threshold_rate import calc_thresholds_rate


import logging
logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def calc_modulation_gain(
        model,
        fms=None,
        cf=10e3,
        model_pars=None,
        m=1,
):
    """Calculate modulation gain of an inner ear model.

    Parameters
    ----------
    fms : array_like
        List of modulation frequencies [Hz].

    TODO: document parameters

    """
    if model_pars is None:
        model_pars = {}

    if fms is None:
        fms = np.logspace(np.log10(10), np.log10(2e3), 16)


    ### Calculate dbspl = threshold + 10 dB
    threshod = calc_thresholds_rate(
        model=model,
        cfs=[cf],
        model_pars=model_pars
    )

    dbspl = threshod.iloc[0] + 10
    log.info("Sound level: {} dB SPL".format(dbspl))

    space = [
        {
            'model': model,
            'fm': fm,
            'cf': cf,
            'dbspl': dbspl,
            'model_pars': model_pars,
            'm': m,
        }
        for fm in fms
    ]

    gains = mr.map(
        _run_model,
        space
    )

    gains = pd.Series(gains, index=fms)
    gains.index.name = 'fm'

    return gains




def _run_model(model, fm, cf, dbspl, model_pars, m):

    duration = 0.6
    onset = 10e-3

    fs = model_pars.setdefault('fs', 100e3)
    model_pars.setdefault('anf_num', (250,0,0))
    model_pars.setdefault('seed', 0)


    sound = wv.am_tone(
        fs=fs,
        fm=fm,
        fc=cf,
        m=m,
        duration=duration,
        dbspl=dbspl,
    )

    anf = model(
        sound=sound,
        cf=cf,
        **model_pars
    )

    si = th.synchronization_index(anf, fm)
    gain = 20 * np.log10(2*si / m)

    return gain
