"""Copyright 2009-2014 Marek Rudnicki

This file is part of cochlea.

cochlea is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

cochlea is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with cochlea.  If not, see <http://www.gnu.org/licenses/>.

"""

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"


import numpy as np
import scipy.signal as dsp

from cochlea.holmberg2007.bm_pars import (
    real_freq_map,
    S_ST,
    S_ED,
    C_eardrum,
    outer_ear_a_48kHz,
    outer_ear_b_48kHz,
    outer_ear_a_100kHz,
    outer_ear_b_100kHz,
    delay_time,
)
from . _traveling_waves import (
    run_bm_wave,
    run_lcr4,
    run_ihcrp,
    run_ihc_meddis2000,
    run_an_sg_carney_holmberg2007
)


# Input signal should be multiplied by this factor
scaling_factor = S_ST * S_ED



def run_middle_ear_filter_orig(signal, fs):
    """Middle ear filter designed using digital wave techinique."""
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


def run_middle_ear_filter(signal, fs):
    """Run middle ear filter model."""

    b,a = _calc_middle_ear_coefs(fs)

    return dsp.lfilter(b, a, signal)


def _calc_outer_ear_coefs(fs):
    if fs == 48000:
        a = outer_ear_a_48kHz
        b = outer_ear_b_48kHz
    elif fs == 100000:
        a = outer_ear_a_100kHz
        b = outer_ear_b_100kHz
    else:
        assert False, "Invalid sampling frequency: fs"

    return b, a


def run_outer_ear_filter(signal, fs):
    """Run outer ear filter."""

    b, a = _calc_outer_ear_coefs(fs)

    return dsp.lfilter(b, a, signal)



def get_nearest_cf(cf):
    """Return nearest valid frequency relative to `cf`."""

    idx = np.argmin(np.abs(real_freq_map - cf))
    nearest_cf = real_freq_map[idx]

    return nearest_cf
