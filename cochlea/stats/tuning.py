#!/usr/bin/env python

"""Cochlear filter tuning.

"""

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd

import mrlib as mr

from . threshold_rate import calc_spont_threshold, calc_threshold


def calc_tuning(
        model,
        cf,
        freqs=None,
        model_pars=None,
        map_backend='serial',
):
    """Calculate runing of the cochlea at `cf`.

    """

    if freqs is None:
        freqs = np.logspace(np.log10(cf/2), np.log10(cf*2), 32)


    spont_rate = mr.apply(
        calc_spont_threshold,
        model=model,
        cf=cf,
        model_pars=model_pars
    )


    space = [
        {
            'model': model,
            'model_pars': model_pars,
            'spont_rate': spont_rate,
            'cf': cf,
            'freq': freq,
        }
        for freq in freqs
    ]


    thresholds = mr.map(
        calc_threshold,
        space,
        backend=map_backend,
    )

    thresholds = pd.Series(thresholds, index=freqs)
    thresholds.index.name = 'freq'


    return thresholds
