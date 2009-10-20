from __future__ import division

import numpy as np
import _catmodel

def run_ihc(sound, cf, fs, cohc=1, cihc=1):
    # TODO: test input parameters; refere to mex files

    # uPa -> Pa
    sound = sound * 1e-6

    vihc = _catmodel.run_ihc(sound, cf, fs, cohc, cihc)

    return vihc



def run_synapse(vihc, cf, fs, anf_type='hsr', powerlaw_implnt='actual', anf_num=1):
    # TODO: test input pars; refere to mex files

    anf_map = {'hsr': 3,
               'msr': 2,
               'lsr': 1}

    implnt_map = {'actual': 1,
                  'approx': 0}

    spike_signal = np.zeros_like( vihc )
    for anf_idx in range(anf_num):
        psth = _catmodel.run_synapse(vihc, cf, fs,
                                     anf_map[anf_type], implnt_map[powerlaw_implnt]);
        spike_signal = spike_signal + psth

    return spike_signal

def set_dB_SPL(dB, signal):
    p0 = 2e-5                   # Pa
    squared = signal**2
    rms = np.sqrt( np.sum(squared) / len(signal) )

    if rms == 0:
        r = 0
    else:
        r = 10**(dB / 20.0) * p0 / rms;

    return signal * r * 1e6     # uPa

