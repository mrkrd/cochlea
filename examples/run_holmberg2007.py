#!/usr/bin/env python

"""Run innear ear model by [Holmberg2007]_.

"""
from __future__ import division, print_function, absolute_import

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
