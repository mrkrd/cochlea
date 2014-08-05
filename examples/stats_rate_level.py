#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2009-2014 Marek Rudnicki

# This file is part of cochlea.

# cochlea is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# cochlea is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with cochlea.  If not, see <http://www.gnu.org/licenses/>.


"""This demo shows show to calculate rate-level characteristic of
a model.

"""
from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"


import numpy as np
import matplotlib.pyplot as plt

import cochlea
from cochlea.stats import calc_rate_level

# import cochlea.external

def main():

    # rates = calc_rate_level(
    #     model=cochlea.run_holmberg2007,
    #     cf=1000,
    #     model_pars={'approximate_cfs': True, 'fs': 48e3}
    #     # dbspls=[0, 20, 50]
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
