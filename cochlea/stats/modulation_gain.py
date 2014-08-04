# Copyright 2009-2014 Marek Rudnicki

# This file is part of cochlea.

# cochlea is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# cochlea is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with cochlea.  If not, see <http://www.gnu.org/licenses/>.

"""Modulation gain of inner ear models.

"""
from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"


import numpy as np
import pandas as pd

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
        map_backend=None,
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
    threshold = calc_thresholds_rate(
        model=model,
        cfs=[cf],
        model_pars=model_pars
    )

    dbspl = threshold['threshold'].iloc[0] + 10
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
    model_pars.setdefault('anf_num', (250,0,0))
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

    si = th.vector_strength(anf, fm)
    gain = 20 * np.log10(2*si / m)

    return gain
