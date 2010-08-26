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



def run_middle_ear_filter(fs, signal):
    """ Middle ear filter model. """
    # Digital wave filter coefficients
    R2_ME = 1. / (2. * fs * C_eardrum)
    R1 = 1. / (2. * np.pi * C_eardrum * 1e3)
    g1_ME = (R2_ME - R1) / (R2_ME + R1)
    Z_ME=0

    # Standard filter coefficients
    b = [(-1-g1_ME), (1+g1_ME)]
    a = [R2_ME, R2_ME*g1_ME]

    return dsp.lfilter(b, a, signal)


def run_outer_ear_filter(fs, signal):
    assert fs == 48000

    a = bm_pars.outer_ear_a
    b = bm_pars.outer_ear_b

    return dsp.lfilter(b, a, signal)


def set_dbspl(dB, signal):
    p0 = 2e-5                   # Pa
    squared = signal**2
    rms = np.sqrt( np.sum(squared) / len(signal) )

    if rms == 0:
        r = 0
    else:
        r = 10**(dB / 20.0) * p0 / rms;

    return signal * r * 1e6     # uPa


