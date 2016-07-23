#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This demo shows how to calculate rate-level characteristic of a
model.

"""
from __future__ import division, absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_rate_level

# import cochlea.external

def main():

    # rates = calc_rate_level(
    #     model=cochlea.run_holmberg2007,
    #     cf=cochlea.get_nearest_cf_holmberg2007(1000),
    #     model_pars={'fs': 48e3, 'syn_mode': 'quantal'},
    # )

    rates = calc_rate_level(
        model=cochlea.run_zilany2014,
        cf=1000,
        model_pars={'fs': 100e3, 'species': 'human'}
    )


    # rates = calc_rate_level(
    #     model=cochlea.external.run_matlab_auditory_periphery,
    #     cf=1000,
    #     model_pars={'fs': 48e3}
    # )

    print(rates)

    rates.plot()

    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
