# Author: Marek Rudnicki
# Time-stamp: <2009-10-06 13:24:05 marek>
#
# Description: Some simple internal tests

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import _pycat

def main():

    fs = 100000.0
    stimdb = 50
    cf = 1000
    cohc = 1
    cihc = 1
    nrep = 100
    t = np.arange(0, 0.05, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = np.sqrt(2) * 20e-6 * 10**(stimdb/20) * np.sin(2*np.pi*cf*t)

    s[0:np.floor(len(s)/2)] = 0

    vihc = _pycat.ihc(s, cf, fs, cohc, cihc)

    plt.plot(s)
    plt.plot(vihc)
    # plt.show()

    synout, psth = _pycat.synapse(vihc, cf, nrep, fs, 3, 0);

    plt.plot(synout)
    plt.show()


def set_dB_SPL(dB, signal):
    p0 = 2e-5                   # Pa
    squared = signal**2
    rms = np.sqrt( np.sum(squared) / len(signal) )

    if rms == 0:
        r = 0
    else:
        r = 10**(dB / 20.0) * p0 / rms;

    return signal * r * 1e6     # uPa



if __name__ == "__main__":
    main()
