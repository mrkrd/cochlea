#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This demo calculates tuning of cochelar filters.

"""

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_tuning


def main():

    ths_human = calc_tuning(
        model=cochlea.run_zilany2013,
        cf=5e3,
        model_pars={'species': 'human'}
    )

    ths_human_glasberg1990 = calc_tuning(
        model=cochlea.run_zilany2013,
        cf=5e3,
        model_pars={'species': 'human_glasberg1990'}
    )

    ths_human.plot(logx=True, label="Shera et al. (2002)")
    ths_human_glasberg1990.plot(logx=True, label="Glasberg & Moore (1990)")

    plt.title("Zilany et al. (2014)")
    plt.legend(loc='best')
    plt.show()


if __name__ == "__main__":
    main()
