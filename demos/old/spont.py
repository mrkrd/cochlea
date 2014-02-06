#!/usr/bin/env python

# Author: Marek Rudnicki
# Time-stamp: <2010-02-18 02:04:48 marek>

# Description:

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

import thorns as th
import pycat

def main():
    fs = 100000                 # Hz
    tmax = 60                   # s
    s = np.zeros(np.ceil(fs*tmax))

    ear = pycat.Zilany2009((1,0,0), powerlaw_implnt='approx')

    anf = ear.run(fs, s)

    folded = th.fold(anf.spikes, 1000)

    rates = [len(train) for train in folded]
    print rates

    plt.plot(rates)
    plt.show()

if __name__ == "__main__":
    main()
