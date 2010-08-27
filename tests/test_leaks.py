#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import thorns.waves as wv
import cochlea

def main():
    fs = 48000
    # fs = 100e3
    s = wv.generate_ramped_tone(fs=fs,
                                freq=1000,
                                tone_duration=10,
                                ramp_duration=1,
                                pad_duration=0,
                                dbspl=60)


    i = 0
    ear = cochlea.Holmberg2007((100,10,10), accumulate=True)
    while True:
        print i
        i = i+1
        # ear = cochlea.Sumner2003((100,10,10), cf=(50,15000,100), accumulate=True)
        anf = ear.run(fs, s)

if __name__ == "__main__":
    main()
