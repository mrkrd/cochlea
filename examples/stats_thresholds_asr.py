#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright 2009-2014 Marek Rudnicki

This file is part of cochlea.

cochlea is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

cochlea is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with cochlea.  If not, see <http://www.gnu.org/licenses/>.


Description
-----------

This demo calculates rate based hearing threshold of an innear ear
model with automatic speach prefiltering in order to adjust thresholds
to human hearing thresholds.

On a fast computer it will take around 1 min to finish.  If you set
shell variable THmap=m, then it will use multiple cores.

"""
from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"


import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_thresholds_rate, calc_human_hearing_thresholds


def main():

    cfs = np.logspace(np.log10(125), np.log10(16e3), 32)

    ths = calc_thresholds_rate(
        model=cochlea.run_zilany2014,
        cfs=cfs,
        model_pars={'species': 'human'},
        asr_filter=True
    )



    human_ths = calc_human_hearing_thresholds(cfs)

    human_ths.plot(logx=True, style='--', linewidth=5)
    ths.plot(logx=True)

    plt.show()


if __name__ == "__main__":
    main()
