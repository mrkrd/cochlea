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

    ### Make sound
    t = np.arange(0, 0.1, 1/fs)
    s = dsp.chirp(t, 80, t[-1], 20000)
    s = cochlea.set_dbspl(s, 50)
    s = np.concatenate( (s, np.zeros(10e-3 * fs)) )



    ### Run model
    anf_trains = cochlea.run_holmberg2007(
        s,
        fs,
        anf_num=(100,0,0),
        seed=0,
    )


    ### Plot auditory nerve response
    anf_acc = th.accumulate(anf_trains, keep=['cf', 'duration'])
    anf_acc.sort('cf', ascending=False, inplace=True)


    fig, ax = plt.subplots(2,1)
    th.plot_signal(
        signal=s,
        fs=fs,
        ax=ax[0]
    )
    th.plot_neurogram(
        anf_acc,
        fs,
        ax=ax[1]
    )
    plt.show()



if __name__ == "__main__":
    main()
