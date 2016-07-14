#!/usr/bin/env python
"""Run inner ear model from [Zilany2014]_.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as dsp

import thorns as th
import cochlea

def main():

    fs = 100e3

    # Make sound
    t = np.arange(0, 0.1, 1/fs)
    s = dsp.chirp(t, 80, t[-1], 20000)
    s = cochlea.set_dbspl(s, 50)
    pad = np.zeros(int(10e-3 * fs))
    sound = np.concatenate( (s, pad) )

    # Run model
    anf = cochlea.run_zilany2014(
        sound,
        fs,
        anf_num=(100,0,0),
        cf=(125, 20000, 100),
        seed=0,
        powerlaw='approximate',
        species='human',
    )

    # Accumulate spike trains
    anf_acc = th.accumulate(anf, keep=['cf', 'duration'])
    anf_acc.sort_values('cf', ascending=False, inplace=True)

    # Plot auditory nerve response
    fig, ax = plt.subplots(2,1)
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
