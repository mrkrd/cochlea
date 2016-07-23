#!/usr/bin/env python

"""Run the external MAP model.

"""
from __future__ import division, absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as dsp

import thorns as th

import cochlea
from cochlea.external import run_matlab_auditory_periphery

def main():

    fs = 48e3

    ### Make sound
    tmax = 0.03
    t = np.arange(0, tmax, 1/fs)
    s = dsp.chirp(t, 80, t[-1], 16000)
    sound = cochlea.set_dbspl(s, 50)



    ### Run model
    anf = run_matlab_auditory_periphery(
        sound,
        fs,
        anf_num=(100,50,20),
        cf=(125, 16000, 80),
        seed=0,
    )



    ### Accumulate spike trains
    anf_acc = th.accumulate(anf, keep=['cf', 'duration'])
    anf_acc.sort('cf', ascending=False, inplace=True)



    ### Plot auditory nerve response
    fig, ax = plt.subplots(2, 1, sharex=True)
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
