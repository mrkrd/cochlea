#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This demo shows show to calculate synchronization of an auditory
model at different characteristic frequencies.

"""
from __future__ import division, absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_synchronization


def main():

    # vss = calc_synchronization(
    #     model=cochlea.run_holmberg2007,
    #     cfs=cochlea.holmberg2007.real_freq_map[10::5],
    #     model_pars={'fs': 48e3}
    # )

    vss = calc_synchronization(
        model=cochlea.run_zilany2014,
        model_pars={'species': 'human'}
    )

    print(vss)

    vss.plot(logx=True)

    plt.show()

if __name__ == "__main__":
    main()
