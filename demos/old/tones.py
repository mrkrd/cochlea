#!/usr/bin/env python

# Author: Marek Rudnicki
# Time-stamp: <2010-02-17 20:45:37 marek>

# Description:

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

import thorns as th
import thorns.waves as wv

import pycat

def main():
    fs = 100000
    cf = 2000
    s = wv.generate_ramped_tone(fs, cf,
                                tone_duration=100,
                                ramp_duration=2.5,
                                pad_duration=100,
                                dbspl=40)

    s = np.tile(s, 200)

    # plt.plot(wv.t(fs,s), s)
    # plt.show()

    print("Generating ANF spikes...")
    ear = pycat.Zilany2009((1,0,0), freq=cf,
                           powerlaw_implnt='approx')
    anf = ear.run(fs, s)
    print("done")

    trains = th.fold(anf.spikes, 200)

    th.plot_raster(trains)
    th.plot_psth(trains, bin_size=1)
    th.plot_isih(trains, bin_size=1)


if __name__ == "__main__":
    main()
