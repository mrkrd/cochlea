# Author: Marek Rudnicki
# Time-stamp: <2009-09-28 20:25:14 marek>
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
    nrep = 1
    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = np.sqrt(2) * 20e-6 * 10**(stimdb/20) * np.sin(2*np.pi*cf*t)

    vich = _pycat.ihc(s, cf, fs, cohc, cihc)

    plt.plot(s)
    plt.plot(vich)
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
