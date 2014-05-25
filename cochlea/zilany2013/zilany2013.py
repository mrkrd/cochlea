#!/usr/bin/env python

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"

import warnings
import itertools
import numpy as np
import pandas as pd

from . import _zilany2013


def run_zilany2013(
        sound,
        fs,
        anf_num,
        cf,
        species,
        seed,
        cohc=1,
        cihc=1,
        powerlaw='approximate',
):
    """Run Zilany et al. (2013/2014) JASA inner ear model.

    Parameters
    ----------
    sound : array_like
        The input sound in Pa.
    fs : float
        Sampling frequency of the sound in Hz.
    anf_num : tuple
        The desired number of auditory nerve fibers per frequency
        channel (CF), (HSR#, MSR#, LSR#).  For example, (100, 75, 25)
        means that we want 100 HSR fibers, 75 MSR fibers and 25 LSR
        fibers per CF.
    cf : float or array_like or tuple
        The center frequency(s) of the simulated auditory nerve fibers.
        If float, then defines a single frequency channel.  If
        array_like (e.g. list or ndarray), then the frequencies are
        used.  If tuple, then must have exactly 3 elements (min_cf,
        max_cf, num_cf) and the frequencies are calculated using the
        Greenwood function.
    species : {'cat', 'human', 'human_glasberg1990'}
        Species.
    seed : int
        Random seed for the spike generator.
    cohc : float <0-1>, optional
        Degredation of the outer hair cells.
    cihc : float <0-1>, optional
        Degredation of the inner hair cells.
    powerlaw : {'approximate', 'actual'}, optional
        Defines which power law implementation should be used.  Can be
        either 'approximate' or 'actual'.

    Returns
    -------
    spike_trains
        Auditory nerve spike trains.

    """
    assert np.max(sound) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1
    assert species in ('cat', 'human', 'human_glasberg1990')

    np.random.seed(seed)

    cfs = _calc_cfs(cf, species)

    channel_args = [
        {
            'signal': sound,
            'cf': cf,
            'fs': fs,
            'cohc': cohc,
            'cihc': cihc,
            'anf_num': anf_num,
            'powerlaw': powerlaw,
            'seed': seed,
            'species': species
        }
        for cf in cfs
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
    seed = args['seed']
    anf_num = args['anf_num']
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
    anf_types = np.repeat(['hsr', 'msr', 'lsr'], anf_num)

    synout = {}
    trains = []
    for anf_type in anf_types:

        if anf_type not in synout:
            ### Run synapse
            synout[anf_type] = _zilany2013.run_synapse(
                fs=fs,
                vihc=vihc,
                cf=cf,
                anf_type=anf_type,
                powerlaw=powerlaw,
                ffGn=False
            )

        ### Run spike generator
        spikes = _zilany2013.run_spike_generator(
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





def _calc_cfs(cf, species):

    if np.isscalar(cf):
        cfs = [float(cf)]

    elif isinstance(cf, tuple) and species == 'cat':
        # Based on GenerateGreenwood_CFList() from DSAM
        # Liberman (1982)
        aA = 456
        k = 0.8
        a = 2.1

        freq_min, freq_max, freq_num = cf

        xmin = np.log10(freq_min / aA + k) / a
        xmax = np.log10(freq_max / aA + k) / a

        x_map = np.linspace(xmin, xmax, freq_num)
        cfs = aA * ( 10**( a*x_map ) - k)

    elif isinstance(cf, tuple) and ('human' in species):
        # Based on GenerateGreenwood_CFList() from DSAM
        # Liberman (1982)
        aA = 165.4
        k = 0.88
        a = 2.1

        freq_min, freq_max, freq_num = cf

        xmin = np.log10(freq_min / aA + k) / a
        xmax = np.log10(freq_max / aA + k) / a

        x_map = np.linspace(xmin, xmax, freq_num)
        cfs = aA * ( 10**( a*x_map ) - k)

    elif isinstance(cf, list) or isinstance(cf, np.ndarray):
        cfs = cf

    else:
        raise RuntimeError("CF must be a scalar, a tuple or a list.")

    return cfs
