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


"""Rate-level charactersitics for inner ear models.

"""
from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"


import numpy as np
import pandas as pd

import thorns as th
import thorns.waves as wv


def calc_rate_level(
        model,
        dbspls=None,
        cf=1000,
        model_pars=None
):
    """Calculate rate-level characteristic of an auditory model.

    """
    if model_pars is None:
        model_pars = {}

    if dbspls is None:
        dbspls = np.arange(-10, 100, 5)

    space = {
        'dbspl': dbspls
    }

    kwargs = {
        'model': model,
        'cf': cf,
        'model_pars': model_pars,
    }

    rates = th.util.map(
        _run_model,
        space,
        kwargs=kwargs,
    )

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
    rate_hsr = th.firing_rate(hsr)

    msr = anf[anf['type']=='msr']
    msr = th.trim(msr, onset, None)
    rate_msr = th.firing_rate(msr)

    lsr = anf[anf['type']=='lsr']
    lsr = th.trim(lsr, onset, None)
    rate_lsr = th.firing_rate(lsr)

    rates = {
        'hsr': rate_hsr,
        'msr': rate_msr,
        'lsr': rate_lsr,
    }

    return rates
