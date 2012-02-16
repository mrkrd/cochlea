#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import thorns.waves as wv
import cochlea

def main():
    fs = 48000
    s = wv.generate_ramped_tone(fs=fs,
                                freq=1000,
                                tone_duration=10,
                                ramp_duration=1,
                                pad_duration=0,
                                dbspl=60)
    for i in range(23):

        print i
        ear = cochlea.Holmberg2007((10,0,10), accumulate=True)
        anf = ear.run(fs, s)

if __name__ == "__main__":
    main()
