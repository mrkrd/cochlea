#!/usr/bin/env python

"""This example illustrates how to run Zilany et al. (2014) model with
the rate (not spikes) as output.

Note that this functionality is given for convenience and not very
well supported.

"""
from __future__ import division, absolute_import, print_function

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
