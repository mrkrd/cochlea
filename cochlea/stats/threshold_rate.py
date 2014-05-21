#!/usr/bin/env python

"""Rate based threshold for inner ear models.

"""

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd
import logging

import mrlib as mr
import mrlib.thorns as th
import mrlib.waves as wv

from . bisection import find_zero

from cochlea.asr import adjust_to_human_thresholds


def calc_thresholds_rate(
        model,
        cfs=None,
        model_pars=None,
        asr_filter=False,
        map_backend='serial'
):
    """Calculate rate based hearing threshold of an inner ear model.

    """
    if cfs is None:
        cfs = np.logspace(np.log10(125), np.log10(16000), 32)

    if model_pars is None:
        model_pars = {}


    ### Calculate spontaneous rate for reference
    spont_rate = mr.apply(
        calc_spont_threshold,
        model=model,
        cf=cfs[0],
        model_pars=model_pars
    )

    space = [
        {
            'model': model,
            'model_pars': model_pars,
            'spont_rate': spont_rate,
            'cf': cf,
            'asr_filter': asr_filter
        }
        for cf in cfs
    ]


    thresholds = mr.map(
        calc_threshold,
        space,
        backend=map_backend,
    )

    thresholds = pd.Series(thresholds, index=cfs)
    thresholds.index.name = 'cf'

    return thresholds




def calc_spont_threshold(model, cf, model_pars):

    pars = dict(model_pars)

    fs = pars.setdefault('fs', 100e3)
    pars.setdefault('seed', 0)
    pars.setdefault('anf_num', (1000,0,0))

    tmax = 250e-3
    silence = np.zeros(fs*tmax)

    anf = model(
        sound=silence,
        cf=cf,
        **pars
    )

    rates = anf.apply(th.rate, axis=1)

    threshold = rates.mean() + rates.std()

    logging.debug("mean + std rate = {}".format(threshold))

    return threshold




def calc_threshold(
        model,
        cf,
        spont_rate,
        model_pars,
        asr_filter=False,
        freq=None
):

    kwargs = {
        'model':model,
        'model_pars': model_pars,
        'cf': cf,
        'spont_rate': spont_rate,
        'asr_filter': asr_filter,
        'freq': freq,
    }

    dbspl_opt = find_zero(
        error_func,
        kwargs=kwargs,
        x1=-100,
        x2=100,
        xtol=0.1
    )

    return dbspl_opt




def error_func(
        dbspl,
        model,
        cf,
        spont_rate,
        model_pars,
        asr_filter=False,
        freq=None
):

    pars = dict(model_pars)

    fs = pars.setdefault('fs', 100e3)
    pars.setdefault('seed', 0)
    pars.setdefault('anf_num', (1000, 0, 0))

    if freq is None:
        freq = cf

    tone_duration = 250e-3
    onset = 15e-3

    sound = wv.ramped_tone(
        fs,
        freq,
        duration=tone_duration,
        ramp=2.5e-3,
        pad=0,
        dbspl=dbspl
    )

    if asr_filter:
        sound = adjust_to_human_thresholds(sound, fs, model)

    anf = model(
        sound=sound,
        cf=cf,
        **pars
    )

    trains = th.trim(anf, onset, None)
    rate = th.rate(trains)

    error = rate - spont_rate

    logging.debug("{} {} {}".format(dbspl, rate, error))

    return error
