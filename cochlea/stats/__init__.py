#!/usr/bin/env python

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"

import numpy as np

from . threshold_rate import calc_thresholds_rate
from . rate_intensity import calc_rate_intensity
from . synchronization import calc_synchronization
from . modulation_gain import calc_modulation_gain


def human_hearing_threshold(freqs):
    """Calculate human hearing thresholds for given frequencies.

    References
    ----------
    E. Terhardt, "Calculating virtual pitch", Hearing Res., vol. 1,
    pp. 155--182, 1979.

    http://www.diracdelta.co.uk/science/source/t/h/threshold%20of%20hearing/source.html

    """
    thresholds = (
        3.64 * freqs**(-0.8)
        - 6.5 * np.exp(-0.6 * (freqs - 3.3)**2)
        + 10**(-3) * freqs**4
    )

    return thresholds
