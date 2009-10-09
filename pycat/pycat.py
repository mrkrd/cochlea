from __future__ import division

import numpy as np
import _pycat

def run_ihc(sound, cf, fs, cohc=1, cihc=1):
    # TODO: test input parameters; refere to mex files

    # uPa -> Pa
    sound = sound * 1e-6

    vihc = _pycat.run_ihc(sound, cf, fs, cohc, cihc)

    return vihc



def run_synapse(vihc, cf, nrep, fs, anf_type='hsr', implnt='actual'):
    # TODO: test input pars; refere to mex files

    anf_map = {'hsr': 3,
               'msr': 2,
               'lsr': 1}

    implnt_map = {'actual': 1,
                  'approx': 0}

    psth = _pycat.run_synapse(vihc, cf, nrep, fs,
                              anf_map[anf_type], implnt_map[implnt]);
    return psth

def set_dB_SPL(dB, signal):
    p0 = 2e-5                   # Pa
    squared = signal**2
    rms = np.sqrt( np.sum(squared) / len(signal) )

    if rms == 0:
        r = 0
    else:
        r = 10**(dB / 20.0) * p0 / rms;

    return signal * r * 1e6     # uPa

