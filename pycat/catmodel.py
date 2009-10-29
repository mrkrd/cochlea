from __future__ import division

import numpy as np
import _catmodel

def run_ihc(sound, cf, fs, cohc=1, cihc=1):
    """
    sound: input sound in uPa
    cf: characteristic frequency
    fs: sampling frequency
    cohc, cihc: degeneration parameters for IHC and OHC cells

    return: IHC receptor potential
    """
    assert isinstance(sound, np.ndarray) and (sound.ndim == 1)
    assert (cf > 80) and (cf < 40e3)
    assert (fs >= 100e3) and (fs <= 500e3)
    assert (cohc >= 0) and (cohc <= 1)
    assert (cihc >= 0) and (cihc <= 1)

    # uPa -> Pa
    # Compatibility with DSAM
    sound = sound * 1e-6

    vihc = _catmodel.run_ihc(sound, cf, fs, cohc, cihc)

    return vihc



def run_synapse(vihc, cf, fs, anf_type='hsr', powerlaw_implnt='actual', anf_num=1):
    """
    vihc: IHC receptor potential
    cf: characteristic frequency
    anf_type: auditory nerve fiber type ('hsr', 'msr' or 'lsr')
    powerlaw_implnt: implementation of the powerlaw ('actual', 'approx')
    anf_num: number of ANF

    return: PSTH from ANF
    """
    assert isinstance(vihc, np.ndarray) and (vihc.ndim == 1)
    assert (cf > 80) and (cf < 40e3)
    assert (fs >= 100e3) and (fs <= 500e3)
    assert anf_type in ['hsr', 'msr', 'lsr']
    assert powerlaw_implnt in ['actual', 'approx']
    assert isinstance(anf_num, int) and (anf_num > 0)

    anf_map = {'hsr': 3,
               'msr': 2,
               'lsr': 1}

    implnt_map = {'actual': 1,
                  'approx': 0}

    spike_signal = np.zeros_like( vihc )
    for anf_idx in range(anf_num):
        psth = _catmodel.run_synapse(vihc, cf, fs,
                                     anf_map[anf_type], implnt_map[powerlaw_implnt])
        spike_signal = spike_signal + psth

    return spike_signal


def set_dB_SPL(dB, signal):
    p0 = 2e-5                   # Pa
    squared = signal**2
    rms = np.sqrt( np.sum(squared) / len(signal) )

    if rms == 0:
        r = 0
    else:
        r = 10**(dB / 20.0) * p0 / rms

    return signal * r * 1e6     # uPa

