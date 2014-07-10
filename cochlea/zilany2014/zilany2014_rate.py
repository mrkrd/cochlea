# Copyright 2013-2014 Marek Rudnicki

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


from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"


import itertools
import numpy as np
import pandas as pd
import itertools

from . import _zilany2014
from . util import calc_cfs

def run_zilany2014_rate(
        sound,
        fs,
        anf_types,
        cf,
        species,
        cohc=1,
        cihc=1,
        powerlaw='approximate',
        ffGn=False
):
    """Run the inner ear model by [Zilany2014]_.  Return mean firing rate
    of the auditory nerve fibers.


    Notes
    -----
    This implementation is was not used very much and may have some
    problems.  Use with caution!  (Like any implementation here, BTW)


    References
    ----------

    .. [Zilany2014] Zilany, M. S., Bruce, I. C., & Carney,
       L. H. (2014). Updated parameters and expanded simulation
       options for a model of the auditory periphery. The Journal of
       the Acoustical Society of America, 135(1), 283-286.

    """
    assert np.max(sound) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1
    assert species in ('cat', 'human')


    if isinstance(anf_types, str):
        anf_types = [anf_types]

    cfs = calc_cfs(cf, species)

    channel_args = [
        {
            'signal': sound,
            'cf': cf,
            'fs': fs,
            'cohc': cohc,
            'cihc': cihc,
            'anf_types': anf_types,
            'powerlaw': powerlaw,
            'species': species,
            'ffGn': ffGn,
        }
        for cf in cfs
    ]


    ### Run model for each channel
    results = map(
        _run_channel,
        channel_args
    )
    results = sum(results, [])


    columns = pd.MultiIndex.from_tuples(
        [(r['anf_type'],r['cf']) for r in results],
        names=['anf_type','cf']
    )
    rates = np.array([r['rate'] for r in results]).T

    rates = pd.DataFrame(
        rates,
        columns=columns
    )

    np.fft.fftpack._fft_cache = {}

    return rates




def _run_channel(args):

    fs = args['fs']
    cf = args['cf']
    signal = args['signal']
    cohc = args['cohc']
    cihc = args['cihc']
    powerlaw = args['powerlaw']
    anf_types = args['anf_types']
    species = args['species']
    ffGn = args['ffGn']


    ### Run BM, IHC
    vihc = _zilany2014.run_ihc(
        signal=signal,
        cf=cf,
        fs=fs,
        species=species,
        cohc=float(cohc),
        cihc=float(cihc)
    )


    duration = len(vihc) / fs


    rates = []
    for anf_type in anf_types:

        ### Run synapse
        synout = _zilany2014.run_synapse(
            fs=fs,
            vihc=vihc,
            cf=cf,
            anf_type=anf_type,
            powerlaw=powerlaw,
            ffGn=ffGn
        )

        rates.append({
            'rate': synout / (1 + 0.75e-3*synout),
            'cf': cf,
            'anf_type': anf_type
        })

    return rates
