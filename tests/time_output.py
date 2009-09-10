# Author: Marek Rudnicki
# Time-stamp: <2009-09-10 14:55:46 marek>
#
# Description: Plot output bitmaps

import numpy as np
import matplotlib.pyplot as plt

import cochlea

def main():
    fs = 48000
    t = np.arange(0, 0.1, 1./fs)
    s = np.sin(2 * np.pi * t * 1000)
    s = cochlea.set_dB_SPL(60, s)

    ear = cochlea.Sumner2002(hsr=100, msr=0, lsr=0,
                             freq=(50, 15000, 100),
                             animal='human')

    hsr, msr, lsr = ear.run(fs, s, times=1, output_format='signals')

    plt.imshow(hsr.T, aspect='auto')
    plt.show()

if __name__ == "__main__":
    main()
