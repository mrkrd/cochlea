#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import thorns as th
import thorns.waves as wv
import thorns.plot as thp

import pycat

def main():

    fs = 100000

    ear = pycat.Zilany2009_Human(anf_num=(100,0,0),
                                 cf=(80, 16000, 100))

    s = wv.generate_ramped_tone(fs,
                                freq=1000,
                                tone_duration=50,
                                ramp_duration=2.5,
                                pad_duration=20,
                                dbspl=70)

    single_trains = ear.run(s, fs)

    anf = th.accumulate(single_trains,
                        ignore=['anf_idx'])

    thp.raster(anf).show()


if __name__ == "__main__":
    main()
