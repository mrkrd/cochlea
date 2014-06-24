#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This demo calculates modulation gain of an innear ear model.

"""

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_modulation_gain


def main():

    gains = calc_modulation_gain(
        model=cochlea.run_zilany2014,
        model_pars={'species': 'cat'}
    )


    # gains = calc_modulation_gain(
    #     model=cochlea.run_holmberg2007,
    #     cf=cochlea.get_nearest_cf_holmberg2007(10e3),
    #     model_pars={'fs': 48e3}
    # )


    # from cochlea.external import run_matlab_auditory_periphery
    # gains = calc_modulation_gain(
    #     model=run_matlab_auditory_periphery,
    #     model_pars={'fs': 48e3}
    # )

    print(gains)

    gains.plot(logx=True)
    plt.show()


if __name__ == "__main__":
    main()
