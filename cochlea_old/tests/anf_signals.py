# Author: Marek Rudnicki
# Time-stamp: <2009-09-15 20:06:40 marek>
#
# Description: Plot output bitmaps

import numpy as np
import matplotlib.pyplot as plt

import cochlea
from traveling_waves import real_freq_map

def main():
    fs = 48000
    t = np.arange(0, 0.1, 1./fs)
    s = np.sin(2 * np.pi * t * 1000)
    s = cochlea.set_dB_SPL(60, s)

    ear = cochlea.Holmberg2008(hsr=100, msr=0, lsr=0)
    hsrH, msr, lsr = ear.run(fs, s, times=1, output_format='signals')



    ear = cochlea.Sumner2002(hsr=100, msr=0, lsr=0,
                             freq=(real_freq_map.min(),
                                   real_freq_map.max(),
                                   100),
                             animal='human')
    hsrS, msr, lsr = ear.run(fs, s, times=1, output_format='signals')



    fig = plt.gcf()
    ax = fig.add_subplot(211)
    ax.set_title("Holmberg 2008")
    ax.imshow(hsrH.T, aspect='auto')
    ax = fig.add_subplot(212)
    ax.set_title("Sumner 2002")
    ax.imshow(hsrS.T, aspect='auto')
    plt.show()

if __name__ == "__main__":
    main()
