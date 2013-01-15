#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import biggles

from mlabwrap import mlab

def main():
    cf = 10000
    cohc = 1
    cihc = 1
    fiberType = 3
    implnt = 0

    f0 = cf
    fs = 100e3
    T = 50e-3
    rt = 5e-3
    stimdb = 10
    nrep = 1

    t = np.arange(0, T, 1/fs)
    s = np.sqrt(2) * 20e-6 * 10**(stimdb/20*np.sin(2*np.pi*f0*t))
    s = np.array([s])
    print s.shape

    vihc = mlab.catmodel_IHC(s, cf, nrep, 1/fs, T*2, cohc, cihc)
    vihc = np.squeeze(vihc)

    print vihc
    biggles.plot(vihc)

if __name__ == "__main__":
    main()
