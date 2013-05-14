#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import warnings
import itertools
import numpy as np
import pandas as pd

import _zilany2013
from zilany2013 import _calc_cfs

def run_zilany2013_rate(
        sound,
        fs,
        anf_types,
        cf,
        species,
        cohc=1,
        cihc=1,
        powerlaw='approximate',
):

    assert np.max(sound) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1
    assert species in ('cat', 'human')


    if isinstance(anf_types, str):
        anf_types = [anf_types]

    cfs = _calc_cfs(cf, species)

    channel_args = [
        {
            'signal': sound,
            'cf': cf,
            'fs': fs,
            'cohc': cohc,
            'cihc': cihc,
            'anf_types': anf_types,
            'powerlaw': powerlaw,
            'species': species
        }
        for cf in cfs
    ]


    ### Run model for each channel
    rates = map(
        _run_channel,
        channel_args
    )


    ### Unpack the results
    rates = pd.DataFrame(rates)
    rates = rates.set_index('cf')

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


    ### Run BM, IHC
    vihc = _zilany2013.run_ihc(
        signal=signal,
        cf=cf,
        fs=fs,
        species=species,
        cohc=float(cohc),
        cihc=float(cihc)
    )


    duration = len(vihc) / fs


    rates = {}
    for anf_type in anf_types:

        ### Run synapse
        synout = _zilany2013.run_synapse(
            fs=fs,
            vihc=vihc,
            cf=cf,
            anf_type=anf_type,
            powerlaw=powerlaw,
            ffGn=False
        )


        rates[anf_type] = synout / (1 + 0.75e-3*synout)


    rates['cf'] = cf

    return rates
