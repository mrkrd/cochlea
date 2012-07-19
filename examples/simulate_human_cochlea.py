#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as dsp

import thorns as th
import thorns.waves as wv

import cochlea

def main():

    fs = 100e3

    ear = cochlea.Zilany2009_Human(
        anf_num=(0,100,0),
        cf=(80, 16000, 100),
        cohc=1.0
    )

    t = np.arange(0, 0.1, 1/fs)
    s = dsp.chirp(t, 80, t[-1], 16000)
    s = cochlea.set_dbspl(s, 50)
    s = np.concatenate( (s, np.zeros(50e-3 * fs)) )

    single_trains = ear.run(s, fs, seed=0)

    anf = th.accumulate(single_trains, ignore=['index'])

    anf_matrix = th.trains_to_signal(anf, fs)

    plt.imshow(anf_matrix.T, aspect='auto')
    plt.show()


if __name__ == "__main__":
    main()
