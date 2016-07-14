#!/usr/bin/env python
"""
Run inner ear model from [Zilany2009]_.

"""
from __future__ import division, absolute_import, print_function


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
    s = np.concatenate((s, np.zeros(int(10e-3*fs))))

    # Run model
    anf = cochlea.run_zilany2009(
        s,
        fs,
        anf_num=(100,0,0),
        cf=(80, 20000, 100),
        seed=0,
        powerlaw='approximate'
    )

    # Plot auditory nerve response
    anf_acc = th.accumulate(anf, keep=['cf', 'duration'])
    anf_acc.sort_values('cf', ascending=False, inplace=True)

    cfs = anf.cf.unique()

    fig, ax = plt.subplots(2,1, sharex=True)
    th.plot_neurogram(
        anf_acc,
        fs,
        ax=ax[0]
    )

    th.plot_raster(
        anf[anf.cf==cfs[30]],
        ax=ax[1]
    )

    ax[1].set_title("CF = {}".format(cfs[30]))

    plt.show()



if __name__ == "__main__":
    main()
