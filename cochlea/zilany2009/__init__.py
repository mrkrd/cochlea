#!/usr/bin/env python

from __future__ import division
from __future__ import print_function


__author__ = "Marek Rudnicki"

import warnings
import itertools
import numpy as np
import pandas as pd

import _pycat


def run_zilany2009(
        sound,
        fs,
        anf_num,
        cf,
        seed,
        cohc=1,
        cihc=1,
        powerlaw='approximate',
):
    """Zilany, M. S. A., Bruce, I. C., Nelson, P. C., and Carney,
    L. H. (2009). A phenomenological model of the synapse between the
    inner hair cell and auditory nerve: long-term adaptation with
    power-law dynamics. The Journal of the Acoustical Society of America,
    126(5):2390-2412.

    """


    assert np.max(sound) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1


    np.random.seed(seed)


    cfs = _calc_cfs(cf)


    ### Run Middle Ear filter
    meout = _pycat.run_middle_ear_filter(
        signal=sound,
        fs=fs
    )


    channel_args = [
        {
            'signal': meout,
            'cf': freq,
            'fs': fs,
            'cohc': cohc,
            'cihc': cihc,
            'anf_num': anf_num,
            'powerlaw': powerlaw,
        }
        for freq in cfs
    ]


    ### Run model for each channel
    nested = map(
        _run_channel,
        channel_args
    )


    ### Unpack the results
    trains = itertools.chain(*nested)
    spike_trains = pd.DataFrame(list(trains))

    np.fft.fftpack._fft_cache = {}

    return spike_trains




def _run_channel(args):

    fs = args['fs']
    cf = args['cf']
    signal = args['signal']
    cohc = args['cohc']
    cihc = args['cihc']
    powerlaw = args['powerlaw']
    anf_num = args['anf_num']


    ### Run BM, IHC
    vihc = _pycat.run_ihc(
        signal=signal,
        cf=cf,
        fs=fs,
        cohc=float(cohc),
        cihc=float(cihc)
    )


    duration = len(vihc) / fs
    anf_types = np.repeat(['hsr', 'msr', 'lsr'], anf_num)

    synout = {}
    trains = []
    for anf_type in anf_types:

        if anf_type not in synout:
            ### Run synapse
            synout[anf_type] = _pycat.run_synapse(
                fs=fs,
                vihc=vihc,
                cf=cf,
                anf_type=anf_type,
                powerlaw=powerlaw,
                ffGn=False
            )

        ### Run spike generator
        spikes = _pycat.run_spike_generator(
            synout=synout[anf_type],
            fs=fs,
        )


        trains.append({
            'spikes': spikes,
            'duration': duration,
            'cf': args['cf'],
            'type': anf_type
        })


    return trains





def _calc_cfs(cf):

    if np.isscalar(cf):
        cfs = [float(cf)]

    elif isinstance(cf, tuple):
        # Based on GenerateGreenwood_CFList() from DSAM
        # Liberman (1982)
        aA = 456
        k = 0.8
        a = 2.1

        freq_min, freq_max, freq_num = cf

        xmin = np.log10( freq_min / aA + k) / a
        xmax = np.log10( freq_max / aA + k) / a

        x_map = np.linspace(xmin, xmax, freq_num)
        cfs = aA * ( 10**( a*x_map ) - k)

    elif isinstance(cf, list) or isinstance(cf, np.ndarray):
        cfs = cf

    else:
        raise RuntimeError("CF must be a scalar, a tuple or a list.")

    return cfs
