#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This demo shows show to calculate rate-intensity characteristic of
a model.

"""

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_rate_intensity


def main():

    # rates = calc_rate_intensity(
    #     model=cochlea.run_holmberg2007,
    #     cf=1000,
    #     model_pars={'approximate_cfs': True, 'fs': 48e3}
    #     # dbspls=[0, 20, 50]
    # )

    rates = calc_rate_intensity(
        model=cochlea.run_zilany2013,
        cf=1000,
        model_pars={'fs': 100e3, 'species': 'human'}
    )

    print(rates)

    rates.plot()

    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
