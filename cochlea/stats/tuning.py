"""
Copyright 2009-2014 Marek Rudnicki

This file is part of cochlea.

cochlea is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

cochlea is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with cochlea.  If not, see <http://www.gnu.org/licenses/>.


Description
-----------
Cochlear filter tuning.

"""
from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"

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
