"""Cochlear filter tuning.

"""

from __future__ import division, print_function, absolute_import

import numpy as np
import pandas as pd

import thorns as th

from . threshold_rate import calc_spont_threshold, calc_threshold


def calc_tuning(
        model,
        cf,
        freqs=None,
        model_pars=None,
        map_backend=None,
):
    """Calculate runing of the cochlea at `cf`."""

    if freqs is None:
        freqs = np.logspace(np.log10(cf/2), np.log10(cf*2), 32)


    spont_rate = th.util.cache(calc_spont_threshold)(
        model=model,
        cf=cf,
        model_pars=model_pars
    )


    space = {
        'freq': freqs
    }

    kwargs = {
        'model': model,
        'model_pars': model_pars,
        'spont_rate': spont_rate,
        'cf': cf,
    }



    thresholds = th.util.map(
        calc_threshold,
        space,
        kwargs=kwargs,
        backend=map_backend,
    )

    thresholds.columns = ['threshold']

    return thresholds
