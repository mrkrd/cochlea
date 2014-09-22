#!/usr/bin/env python

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

"""Run innear ear model by [Holmberg2007]_.

"""
from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as dsp

import thorns as th

import cochlea

def main():

    fs = 48e3
    tmax = 0.1

    ### Make sound
    t = np.arange(0, tmax, 1/fs)
    s = np.zeros_like(t)


    ### Run model
    vesicle_trains = cochlea.run_holmberg2007_vesicles(
        s,
        fs,
        anf_num=(1,0,0),
        seed=0,
    )



    print(vesicle_trains)


    ### Calculate average rate
    all_vesicles = np.concatenate(vesicle_trains['vesicles'])
    rate = len(all_vesicles) / tmax / len(vesicle_trains['vesicles'])
    print()
    print("Rate:", rate)



if __name__ == "__main__":
    main()
