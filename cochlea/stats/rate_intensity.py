#!/usr/bin/env python

"""Rate-intensity charactersitics for inner ear models.

"""

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd

import thorns as th
import thorns.waves as wv


def calc_rate_intensity(
        model,
        dbspls=None,
        cf=1000,
        model_pars=None
):
    """Calculate rate-intensity characteristic of an auditory model.

    """
    if model_pars is None:
        model_pars = {}

    if dbspls is None:
        dbspls = np.arange(-10, 100, 5)

    space = [
        {
            'model': model,
            'dbspl': dbspl,
            'cf': cf,
            'model_pars': model_pars,
        }
        for dbspl in dbspls
    ]

    rates = th.util.map(
        _run_model,
        space
    )

    rates = pd.DataFrame(list(rates))
    rates = rates.set_index('dbspl')

    return rates



def _run_model(model, dbspl, cf, model_pars):

    duration = 100e-3
    onset = 10e-3

    fs = model_pars.setdefault('fs', 100e3)
    model_pars.setdefault('anf_num', (250,250,250))
    model_pars.setdefault('seed', 0)

    sound = wv.ramped_tone(
        fs=fs,
        freq=cf,
        duration=duration,
        pad=0,
        dbspl=dbspl
    )

    anf = model(
        sound=sound,
        cf=cf,
        **model_pars
    )

    hsr = anf[anf['type']=='hsr']
    hsr = th.trim(hsr, onset, None)
    rate_hsr = th.rate(hsr)

    msr = anf[anf['type']=='msr']
    msr = th.trim(msr, onset, None)
    rate_msr = th.rate(msr)

    lsr = anf[anf['type']=='lsr']
    lsr = th.trim(lsr, onset, None)
    rate_lsr = th.rate(lsr)

    rates = {
        'dbspl': dbspl,
        'hsr': rate_hsr,
        'msr': rate_msr,
        'lsr': rate_lsr
    }

    return rates
