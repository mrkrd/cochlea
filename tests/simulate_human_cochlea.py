#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import matplotlib.pyplot as plt

import thorns as th
import thorns.waves as wv

import pycat

def main():

    fs = 100000

    ear = pycat.Zilany2009_Human(anf_num=(100,0,0),
                                 cf=(80, 16000, 100))

    s = wv.generate_ramped_tone(fs,
                                freq=600,
                                tone_duration=50,
                                ramp_duration=2.5,
                                pad_duration=20,
                                dbspl=50)

    single_trains = ear.run(s, fs)

    anf = th.accumulate(single_trains,
                        ignore=['anf_idx'])

    anf_matrix = th.trains_to_signal(anf, fs)
    plt.imshow(anf_matrix.T, aspect='auto')
    plt.show()


if __name__ == "__main__":
    main()
