#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This demo calculates tuning of cochlear filters.

On a fast computer it will take around 2 minutes to finish.  If you
set shell variable THmap=m, then it will use multiple cores.

"""
from __future__ import division, absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_tuning


def main():

    ths_human = calc_tuning(
        model=cochlea.run_zilany2014,
        cf=5e3,
        model_pars={'species': 'human'}
    )

    print(ths_human)


    fig,ax = plt.subplots()

    ths_human.plot(logx=True, ax=ax)

    ax.set_title("Zilany et al. (2014)")
    ax.legend(loc='best')

    plt.show()


if __name__ == "__main__":
    main()
