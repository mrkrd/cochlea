#!/usr/bin/env python

# Author: Marek Rudnicki
# Time-stamp: <2009-12-14 15:35:27 marek>

# Description:

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

import thorns as th
import thorns.waves as wv

import pycat

def main():
    fs = 100000
    cf = 15000
    s = wv.generate_ramped_tone(fs, cf,
                                tone_duration=50,
                                ramp_duration=2.5,
                                pad_duration=55,
                                dbspl=50)

    s = np.tile(s, 50)

    # plt.plot(wv.t(fs,s), s)
    # plt.show()

    print("Generating ANF spikes...")
    ear = pycat.Zilany2009(hsr=1, msr=0, lsr=0, freq=cf,
                             powerlaw_implnt='actual')
    hsr, msr, lsr = ear.run(fs, s, times=1)
    print("done")

    hsr_trains = th.fold( hsr['spikes'], 105 )

    th.plot_raster(hsr_trains)
    th.plot_psth(hsr_trains, bin_size=0.5)
    th.plot_isih(hsr_trains, bin_size=0.5)


if __name__ == "__main__":
    main()
