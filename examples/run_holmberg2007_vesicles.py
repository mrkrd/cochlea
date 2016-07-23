#!/usr/bin/env python

"""Run innear ear model by [Holmberg2007]_ in quantal mode and output
vesicle events instead of spikes.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as dsp

import thorns as th

import cochlea

def main():

    fs = 48e3
    tmax = 0.1

    ### Make sound
    t = np.arange(0, tmax, 1/fs)
    s = np.zeros_like(t)


    ### Run model
    vesicle_trains = cochlea.run_holmberg2007_vesicles(
        s,
        fs,
        anf_num=(1,0,0),
        seed=0,
    )



    print(vesicle_trains)

    ### Need to rename a column: vesicles -> spikes, because that's
    ### the expected name by many functions in thorns
    trains = vesicle_trains.rename(columns={'vesicles': 'spikes'})

    ### Calculate average rate
    rate = th.firing_rate(trains)
    print()
    print("Spontanious rate:", rate)


    ### Raster plot
    th.plot_raster(trains)
    th.show()


if __name__ == "__main__":
    main()
