#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Automatic Speach Recognition related functions.

"""

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd
import scipy

import inspect
import os


def adjust_to_human_thresholds(signal, fs, model):
    """Filter `signal` in order to adjust threshold of `model` to human
    hearing thresholds.

    """

    if inspect.isfunction(model):
        name = model.func_name
    else:
        name = str(model)


    if 'zilany2014' in name:
        fname = _data_dir("zilany2014.csv")
    elif 'holmberg2007' in name:
        fname = _data_dir("holmberg2007.csv")
    elif 'matlab_auditory_periphery' in name:
        fname = _data_dir("matlab_auditory_periphery.csv")
    else:
        raise NotImplementedError


    filt = pd.read_csv(fname)

    signal_fft = np.fft.rfft(signal)
    freqs = np.linspace(0, fs/2, len(signal_fft))


    # Interpolate the filter to fit signal's FFT
    h = scipy.interpolate.interp1d(
        x=filt.freq,
        y=filt.h,
        kind='cubic',
        bounds_error=False,
        fill_value=-100
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
