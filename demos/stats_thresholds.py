#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This demo calculates rate based hearing threshold of an innear ear
model.

"""

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_thresholds_rate


def main():

    ths = calc_thresholds_rate(
        model=cochlea.run_zilany2013,
        model_pars={'species': 'cat'}
    )
    print(ths)

    ths.plot(logx=True)
    plt.show()


if __name__ == "__main__":
    main()
