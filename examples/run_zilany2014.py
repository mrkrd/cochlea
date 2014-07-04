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


"""Run inner ear model from [Zilany2014]_.

"""
from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"


import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as dsp

import thorns as th

import cochlea

def main():

    fs = 100e3

    ### Make sound
    t = np.arange(0, 0.1, 1/fs)
    s = dsp.chirp(t, 80, t[-1], 20000)
    s = cochlea.set_dbspl(s, 50)
    pad = np.zeros(10e-3 * fs)
    sound = np.concatenate( (s, pad) )



    ### Run model
    anf = cochlea.run_zilany2014(
        sound,
        fs,
        anf_num=(100,0,0),
        cf=(125, 20000, 100),
        seed=0,
        powerlaw='approximate',
        species='human',
    )


    ### Accumulate spike trains
    anf_acc = th.accumulate(anf, keep=['cf', 'duration'])
    anf_acc.sort('cf', ascending=False, inplace=True)



    ### Plot auditory nerve response
    fig, ax = plt.subplots(2,1)
    th.plot_signal(
        signal=sound,
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
