#!/usr/bin/env python
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

This example illustrates how to run Zilany et al. (2014) model with
the rate (not spikes) as output.

Note that this functionality is given for convenience and not very
well supported.

"""
from __future__ import division, absolute_import, print_function

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
    s = np.concatenate( (s, np.zeros(10e-3 * fs)) )



    ### Run model
    rates = cochlea.run_zilany2014_rate(
        s,
        fs,
        anf_types=['msr'],
        cf=(125, 20000, 100),
        powerlaw='approximate',
        species='human'
    )


    ### Plot rates
    fig, ax = plt.subplots()
    img = ax.imshow(
        rates.T,
        aspect='auto'
    )
    plt.colorbar(img)
    plt.show()



if __name__ == "__main__":
    main()
