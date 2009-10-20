# Author: Marek Rudnicki
# Time-stamp: <2009-10-21 01:21:25 marek>
#
# Description: Plot output bitmaps

import numpy as np
import matplotlib.pyplot as plt

import pycat

def main():
    fs = 100000
    t = np.arange(0, 0.5, 1./fs)
    s = np.sin(2 * np.pi * t * 1000)
    s = pycat.set_dB_SPL(60, s)

    ear = pycat.Carney2009(hsr=50, msr=0, lsr=0,
                           freq=(50, 10000, 100),
                           powerlaw_implnt='approx')
    hsr, msr, lsr = ear.run(fs, s, times=1, output_format='signals')

    ax = plt.gca()
    ax.imshow(hsr.T, aspect='auto')
    plt.show()

if __name__ == "__main__":
    main()
