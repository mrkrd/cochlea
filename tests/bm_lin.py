# Author: Marek Rudnicki
# Time-stamp: <2009-09-17 12:51:25 marek>
#
# Description: Testing input / output of BM

import numpy as np
import matplotlib.pyplot as plt

import dsam
import traveling_waves as tw
from cochlea.auditory_periphery import par_dir


def main():
    fs = 48000.

    s = np.zeros( np.floor(0.1 * fs) )
    idx = np.floor( len(s) / 4 )
    s[idx] = 100

    bm = dsam.EarModule("BM_DRNL")
    bm.read_pars(par_dir("bm_drnl_gp.par"))
    bm.set_par("CF_MODE", "single")
    bm.set_par("SINGLE_CF", 1000)

    bm.run(fs, s)
    ref = bm.get_signal()

    bm.run(fs, s * 2)
    double = bm.get_signal()

    fig = plt.gcf()

    ax = fig.add_subplot(211)
    ax.plot(ref)
    ax.plot(double / 2.)


    tw_ref = tw.run_bm(fs, s)[:,60]
    tw_double = tw.run_bm(fs, s*2)[:,60]
    ax2 = fig.add_subplot(212)
    ax2.plot(tw_ref)
    ax2.plot(tw_double / 2)

    plt.show()

if __name__ == "__main__":
    main()
