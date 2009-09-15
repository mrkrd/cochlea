# Author: Marek Rudnicki
# Time-stamp: <2009-09-14 19:31:42 marek>
#
# Description: Computes traveling wave delays (signal-front delays) of
# basilar membrane models.

import numpy as np
import matplotlib.pyplot as plt

import cochlea

def delay_vs_intensity_holmberg2008():
    ear = cochlea.Holmberg2008(hsr=1, msr=0, lsr=0)

    fs = 48000.0
    click = np.zeros( np.ceil(0.1 * fs) )
    idx = np.floor( len(click) / 4)

    amp_list = [1, 10]

    for amp in amp_list:
        click[idx] = amp

        ear.run(fs, click)

        bm_delay = np.argmax(np.abs(ear.bm_signal) > 1e-12, axis=0)

        plt.plot(bm_delay - idx)

    plt.show()


def delay_vs_intensity_sumner2002():
    ear = cochlea.Sumner2002(hsr=1, msr=0, lsr=0)

    fs = 100000
    click = np.zeros( np.ceil(0.1 * fs) )
    idx = np.floor( len(click) / 4)

    amp_list = [1, 10]

    for amp in amp_list:
        click[idx] = amp

        ear.run(fs, click)

        bm_delay = np.argmax(np.abs(ear.bm.get_signal) != 0, axis=0)

        plt.plot(bm_delay - idx)

    plt.show()


if __name__ == "__main__":
    delay_vs_intensity_sumner2002()
