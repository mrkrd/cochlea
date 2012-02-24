# Author: Marek Rudnicki
# Time-stamp: <2009-09-17 11:40:13 marek>
#
# Description: IHC RP in response to tone

import numpy as np
import matplotlib.pyplot as plt

import cochlea

def main():
    cf = 10000

    fs = 100000.
    t = np.arange(0, 0.1, 1./fs)
    s = np.sin(2 * np.pi * t * cf)
    s = cochlea.set_dB_SPL(80, s)

    earS = cochlea.Sumner2002(hsr=0, msr=0, lsr=0, freq=cf)
    earL = cochlea.LopezPoveda2006(hsr=0, msr=0, lsr=0, freq=cf)

    earS.run(fs, s)
    earL.run(fs, s)

    rpS = earS.bm.get_signal() - earS.bm.get_signal().mean()
    rpL = earL.bm.get_signal() - earL.bm.get_signal().mean()

    ax = plt.gca()
    ax.plot(rpS)
    ax.plot(rpL)
    plt.show()


if __name__ == "__main__":
    main()
