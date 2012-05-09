#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import matplotlib.pyplot as plt

import thorns as th
import thorns.waves as wv

import cochlea

def main():

    fs = 100e3

    ear = cochlea.Zilany2009_Human(
        anf_num=(0,1000,0),
        cf=(80, 16000, 100),
        cohc=1.0
    )

    s = wv.generate_ramped_tone(
        fs,
        freq=600,
        tone_duration=50e-3,
        ramp_duration=2.5e-3,
        pad_duration=20e-3,
        dbspl=60
    )
    s = np.concatenate( (s, np.zeros(50e-3 * fs)) )

    single_trains = ear.run(s, fs, seed=0)

    anf = th.accumulate(single_trains, ignore=['index'])

    anf_matrix = th.trains_to_signal(anf, fs)

    plt.imshow(anf_matrix.T, aspect='auto')
    plt.show()


if __name__ == "__main__":
    main()
