# Description: Model of auditory periphery of: Zilany, M.S.A., Bruce,
# I.C., Nelson, P.C., and Carney, L.H. (manuscript in preparation) 2009

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import scipy.interpolate
import pandas as pd
import os
import itertools

import _pycat
from zilany2009 import _run_channel


def run_zilany2009_human(
        sound,
        fs,
        anf_num,
        seed,
        cf,
        cohc=1,
        cihc=1,
        powerlaw='approximate',
        middle_ear=True,
):




    assert np.max(sound) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1


    np.random.seed(seed)


    cfs = _calc_cfs(cf)



    ### Run human middle ear filter
    if middle_ear:
        meout = _run_human_me_filter_for_zilany2009(
            signal=sound,
            fs=fs
        )
    else:
        meout = sound


    channel_args = [
        {
            'signal': meout,
            'cf': freq,
            'fs': fs,
            'cohc': cohc,
            'cihc': cihc,
            'anf_num': anf_num,
            'powerlaw': powerlaw,
            'seed': seed
        }
        for freq in cfs
    ]


    ### Run the model for each channel
    nested = map(
        _run_channel,
        channel_args
    )


    ### Unpack the results
    trains = itertools.chain(*nested)
    spike_trains = pd.DataFrame(
        list(trains)
    )

    np.fft.fftpack._fft_cache = {}

    return spike_trains





def _calc_cfs(cf):

    if np.isscalar(cf):
        cfs = [float(cf)]

    elif isinstance(cf, tuple):
        # Based on GenerateGreenwood_CFList() from DSAM
        # Liberman (1982)
        aA = 165.4
        k = 0.88
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









def _run_human_me_filter_for_zilany2009(signal, fs):
    assert fs > 40e3

    signal_fft = np.fft.rfft(signal)
    freqs = np.linspace(0, fs/2, len(signal_fft))

    me_filter = pd.read_csv(
        _data_dir('human_me_filter.csv')
    )



    # Interpolate the filter to fit signal's FFT
    h = scipy.interpolate.interp1d(
        x=me_filter.freq,
        y=me_filter.h,
        kind='cubic',
        bounds_error=False,
        fill_value=0
    )



    # Convert dB to amplitude ratio
    h_ratio = 10 ** (h(freqs) / 20)



    # Apply the filter
    signal_fft *= h_ratio


    # Go back to time domain
    filtered = np.fft.irfft(signal_fft)


    return filtered


def _data_dir(par_file):
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, par_file)
