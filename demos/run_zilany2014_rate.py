#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as dsp

import mrlib.thorns as th

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
