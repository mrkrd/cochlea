# Author: Marek Rudnicki
# Time-stamp: <2009-09-15 19:33:51 marek>
#
# Description: Computes traveling wave delays (signal-front delays) of
# basilar membrane models.

import numpy as np
import matplotlib.pyplot as plt

import cochlea
import traveling_waves as tw

def delay_vs_intensity_holmberg2008():
    ear = cochlea.Holmberg2008(hsr=1, msr=0, lsr=0)

    fs = 48000
    click = np.zeros( np.ceil(0.1 * fs) )
    idx = np.floor( len(click) / 4)

    amp_list = [1]

    for amp in amp_list:
        click[idx] = amp

        ear.run(fs, click)

        bm_delay = np.argmax(np.abs(ear.bm_signal), axis=0)

        plt.plot(bm_delay - idx)



def delay_vs_intensity_sumner2002():

    ear = cochlea.Sumner2002(hsr=1, msr=0, lsr=0,
                             freq=(tw.real_freq_map.min(),
                                   tw.real_freq_map.max(),
                                   100))

    fs = 100000
    click = np.zeros( np.ceil(0.1 * fs) )
    idx = np.floor( len(click) / 4)

    amp_list = [1]

    for amp in amp_list:
        click[idx] = amp

        ear.run(fs, click)

        bm_delay = np.argmax(np.abs(ear.bm.get_signal()), axis=0)

        plt.plot(bm_delay - idx)




if __name__ == "__main__":
    delay_vs_intensity_sumner2002()
    delay_vs_intensity_holmberg2008()
    plt.show()
