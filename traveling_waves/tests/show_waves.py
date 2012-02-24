# Author: Marek Rudnicki
# Time-stamp: <2010-05-22 00:12:43 marek>
#
# Description: Show click response of the BM model


import numpy as np
import matplotlib.pyplot as plt

import traveling_waves as tw

def main():
    fs = 48000.
    s = np.zeros( np.floor(0.1 * fs) )
    idx = np.floor( len(s) / 4. )
    s[idx] = 1

    xbm = tw.run_bm(fs, s, mode='x', enable_LCR4=True)

    plt.imshow(xbm, aspect='auto')
    plt.show()

if __name__ == "__main__":
    main()

