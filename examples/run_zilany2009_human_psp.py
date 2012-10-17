#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as dsp

import cochlea

def main():

    fs = 100e3

    ### Make sound
    t = np.arange(0, 0.1, 1/fs)

    s = dsp.chirp(t, 80, t[-1], 20000)
    s = np.concatenate( (s, np.zeros(10e-3 * fs)) )

    s = cochlea.set_dbspl(s, 50)



    ### Run model
    rates, cfs = cochlea.run_zilany2009_human_psp(
        s,
        fs,
        anf_type='msr',
        cf=(80, 20000, 100),
    )




    ### Plot
    fig, ax = plt.subplots()
    img = ax.imshow(
        rates.T,
        aspect='auto'
    )

    plt.colorbar(img)


    plt.show()



if __name__ == "__main__":
    main()
