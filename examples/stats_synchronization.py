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

This demo shows show to calculate synchronization of an auditory
model at different characteristic frequencies.

"""
from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"


import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_synchronization


def main():

    # sis = calc_synchronization(
    #     model=cochlea.run_holmberg2007,
    #     cfs=cochlea.holmberg2007.real_freq_map[10::5],
    #     model_pars={'fs': 48e3}
    # )

    sis = calc_synchronization(
        model=cochlea.run_zilany2014,
        model_pars={'species': 'human'}
    )

    print(sis)

    hsr_sis = sis.pivot(index='dbspl', columns='cf', values='hsr')

    fig,ax = plt.subplots(2,1)

    ax[0].imshow(hsr_sis, aspect='auto', interpolation='nearest')
    hsr_sis.max().plot(ax=ax[1], logx=True)

    plt.show()

if __name__ == "__main__":
    main()
