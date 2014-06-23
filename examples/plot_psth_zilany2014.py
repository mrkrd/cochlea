#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Run Zilany et al. (2014) model and plot a PSTH.

Tip: This script has no extra dependencies (e.g. thorns).  However,
it's recommended to install thorns and use plot_psth() instead of
using own or even this this PSTH implementation.

"""

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt

import cochlea

def main():
    fs = 100e3
    cf = 8e3

    ### Make sound
    t = np.arange(0, 0.1, 1/fs)
    tone = np.sin(2*np.pi*t*cf)

    tone = cochlea.set_dbspl(tone, 20)

    pad = np.zeros(50e-3 * fs)
    sound = np.concatenate( (tone, pad) )



    ### Run model
    anf = cochlea.run_zilany2014(
        sound,
        fs,
        anf_num=(200,0,0),
        cf=cf,
        seed=0,
        powerlaw='approximate',
        species='human'
    )


    print(anf.head(20))


    ### Plot PSTH
    all_spikes = np.concatenate(anf.spikes)
    tmax = anf.duration.max()

    bin_size = 1e-3

    fig,ax = plt.subplots()

    ax.hist(
        all_spikes,
        bins=int(tmax/bin_size),
        range=(0,tmax),
        weights=np.ones_like(all_spikes)/bin_size/len(anf)
    )

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Rate [spikes/s]")
    ax.set_title("PSTH")

    plt.show()


if __name__ == "__main__":
    main()
