from __future__ import division

import numpy as np
import scipy.signal as dsp

import bm_pars
from bm_pars import real_freq_map, S_ST, S_ED, C_eardrum
from _tw import *


# Input signal should be multiplied by this factor
scaling_factor = S_ST * S_ED



def run_middle_ear_filter_orig(fs, signal):
    """ Middle ear filter designed using digital wave techinique. """
    assert fs == 48000

    R2_ME = 1. / (2. * fs * C_eardrum)
    R1 = 1. / (2. * np.pi * C_eardrum * 1e3)

    g1_ME = (R2_ME - R1) / (R2_ME + R1)
    Z_ME=0

    out = np.zeros(len(signal))

    for i,samp in enumerate(signal):
        b2 = samp + g1_ME * (samp - Z_ME)
        Z_ME_b = Z_ME
        Z_ME = b2
        out[i] = ((Z_ME_b - b2) / R2_ME)

    return out


def _calc_middle_ear_coefs(fs):
    # Digital wave filter coefficients
    R2_ME = 1. / (2. * fs * C_eardrum)
    R1 = 1. / (2. * np.pi * C_eardrum * 1e3)
    g1_ME = (R2_ME - R1) / (R2_ME + R1)
    Z_ME=0

    # Standard filter coefficients
    b = [(-1-g1_ME), (1+g1_ME)]
    a = [R2_ME, R2_ME*g1_ME]

    return b, a


def run_middle_ear_filter(fs, signal):
    """ Middle ear filter model. """
    b,a = _calc_middle_ear_coefs(fs)

    return dsp.lfilter(b, a, signal)


def _calc_outer_ear_coefs(fs):
    if fs == 48000:
        a = bm_pars.outer_ear_a_48kHz
        b = bm_pars.outer_ear_b_48kHz
    elif fs == 100000:
        a = bm_pars.outer_ear_a_100kHz
        b = bm_pars.outer_ear_b_100kHz
    else:
        assert False, "Invalid sampling frequency: fs"

    return b, a


def run_outer_ear_filter(fs, signal):
    b, a = _calc_outer_ear_coefs(fs)

    return dsp.lfilter(b, a, signal)



def find_closest_freq_idx_in_map(freq):
    m = np.abs(real_freq_map - freq)

    return int(np.argmin(m))


def find_closest_freq_in_map(freq):
    idx = find_closest_freq_idx_in_map(freq)
    return real_freq_map[idx]
