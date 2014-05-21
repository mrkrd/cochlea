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

    ths = calc_tuning(
        model=cochlea.run_zilany2013,
        cf=5e3,
        model_pars={'species': 'human'}
    )

    ths.plot(logx=True)
    plt.show()


if __name__ == "__main__":
    main()
