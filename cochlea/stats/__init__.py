#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np

from threshold_rate import calc_hearing_thresholds_rate
from threshold_si import calc_hearing_thresholds_si

def calc_human_hearing_threshold(freqs):
    """Terhardt, E. (1979). Calculating virtual pitch. Hearing
    Research, 1(2):155-182.

    http://www.diracdelta.co.uk/science/source/t/h/threshold%20of%20hearing/source.html

    """
    f = freqs / 1000                # kHz -> Hz

    th = 3.64 * f**(-0.8) - \
         6.5 * np.exp(-0.6 * (f - 3.3)**2) + \
         10**(-3) * f**4

    return th
