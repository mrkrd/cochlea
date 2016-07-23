"""Modulation gain of inner ear models.

"""
from __future__ import division, print_function, absolute_import

import numpy as np

import thorns as th
import thorns.waves as wv

from . threshold_rate import calc_thresholds_rate

import logging

log = logging.getLogger('thorns')


def calc_modulation_gain(
        model,
        fms=None,
        cf=10e3,
        model_pars=None,
        m=1,
        level_above_threshold=10,
        map_backend=None,
):
    """Calculate modulation gain of an inner ear model.

    Parameters
    ----------
    fms : array_like
        List of modulation frequencies (Hz).
    cf : scalar, optional
        Characteristic frequency (Hz).
    model_pars : dict, ooptional
        Additional parameters to be passed to the model funtion.
    m : float, optional
        Modulation depth in the range of <0, 1>.
    level_above_threshold : scalar, optional
        Sound level which will be added to the threshold level to
        calculate modulation gain.

    map_backend : str, optional
        Map backend that will be used in the calculation.  See
        `thorns.util.map()`.

    """
    if model_pars is None:
        model_pars = {}

    if fms is None:
        fms = np.logspace(np.log10(10), np.log10(2e3), 16)

    # Calculate dbspl = threshold + 10 dB
    threshold = calc_thresholds_rate(
        model=model,
        cfs=[cf],
        model_pars=model_pars
    )

    dbspl = threshold['threshold'].iloc[0] + level_above_threshold
    log.info("Sound level: {} dB SPL".format(dbspl))

    space = {
        'fm': fms,
    }

    kwargs = {
        'model': model,
        'cf': cf,
        'dbspl': dbspl,
        'model_pars': model_pars,
        'm': m,
    }

    gains = th.util.map(
        _run_model,
        space,
        kwargs=kwargs,
        backend=map_backend,
    )

    gains.columns = ['gain']

    return gains


def _run_model(model, fm, cf, dbspl, model_pars, m):

    duration = 0.6
    onset = 10e-3

    fs = model_pars.setdefault('fs', 100e3)
    model_pars.setdefault('anf_num', (250, 0, 0))
    model_pars.setdefault('seed', 0)

    sound = wv.amplitude_modulated_tone(
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

    trimmed = th.trim(anf, onset)

    vs = th.vector_strength(trimmed, fm)
    gain = 20 * np.log10(2*vs / m)

    return gain
