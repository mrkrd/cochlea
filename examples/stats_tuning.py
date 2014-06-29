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

This demo calculates tuning of cochelar filters.

On a fast computer it will take around 2 minutes to finish.  If you
set shell variable THmap=m, then it will use multiple cores.

"""
from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"


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

    ths_human_glasberg1990 = calc_tuning(
        model=cochlea.run_zilany2014,
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
