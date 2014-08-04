#!/usr/bin/env python

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd

from . threshold_rate import calc_thresholds_rate
from . rate_level import calc_rate_level
from . synchronization import calc_synchronization
from . modulation_gain import calc_modulation_gain
from . tuning import calc_tuning


def calc_human_hearing_thresholds(freqs):
    """Calculate human hearing thresholds for given frequencies.

    References
    ----------
    E. Terhardt, "Calculating virtual pitch", Hearing Res., vol. 1,
    pp. 155--182, 1979.

    http://www.diracdelta.co.uk/science/source/t/h/threshold%20of%20hearing/source.html

    """
    f = freqs/1000

    thresholds = (
        3.64 * f**(-0.8)
        - 6.5 * np.exp(-0.6 * (f - 3.3)**2)
        + 10**(-3) * f**4
    )

    thresholds = pd.Series(thresholds, index=freqs)

    return thresholds
