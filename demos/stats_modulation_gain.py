#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This demo calculates modulation gaoin of an innear ear model.

"""

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_modulation_gain


def main():

    gains = calc_modulation_gain(
        model=cochlea.run_zilany2013,
        model_pars={'species': 'cat'}
    )

    print(gains)

    gains.plot(logx=True)
    plt.show()


if __name__ == "__main__":
    main()
