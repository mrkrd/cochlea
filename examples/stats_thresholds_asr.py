#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This demo calculates rate based hearing threshold of an innear ear
model with automatic speach prefiltering in order to adjust thresholds
to human hearing thresholds.

On a fast computer it will take around 1 min to finish.  If you set
shell variable THmap=m, then it will use multiple cores.

"""
from __future__ import division, absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt

import cochlea.external
from cochlea.stats import calc_thresholds_rate, calc_human_hearing_thresholds


def main():

    cfs = np.logspace(np.log10(125), np.log10(16e3), 32)

    ths = calc_thresholds_rate(
        model=cochlea.run_zilany2014,
        cfs=cfs,
        model_pars={'species': 'human'},
        asr_filter=True
    )


    # ths = calc_thresholds_rate(
    #     model=cochlea.external.run_matlab_auditory_periphery,
    #     cfs=cfs,
    #     model_pars={},
    #     asr_filter=True
    # )


    human_ths = calc_human_hearing_thresholds(cfs)

    fig, ax = plt.subplots()

    human_ths.plot(ax=ax, logx=True, style='--', linewidth=5)
    ths.plot(ax=ax, logx=True)

    plt.show()


if __name__ == "__main__":
    main()
