#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as dsp

import marlib.thorns as th

import cochlea

def main():

    fs = 100e3

    ### Make sound
    t = np.arange(0, 0.1, 1/fs)
    s = dsp.chirp(t, 80, t[-1], 20000)
    s = cochlea.set_dbspl(s, 50)
    s = np.concatenate( (s, np.zeros(10e-3 * fs)) )



    ### Run model
    anf = cochlea.run_zilany2009(
        s,
        fs,
        anf_num=(100,0,0),
        cf=(80, 20000, 100),
        seed=0,
        parallel=True
    )



    ### Plot auditory nerve response
    anf_acc = th.accumulate(anf, keep=['cf', 'duration'])
    anf_acc.sort('cf', ascending=False, inplace=True)


    fig, ax = plt.subplots()
    th.plot_neurogram(
        anf_acc,
        fs,
        axis=ax
    )
    plt.show()


if __name__ == "__main__":
    main()
