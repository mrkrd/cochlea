# PSTH at the threshold intensity at low frequency

import numpy as np
import matplotlib.pyplot as plt

import thorns as th
import cochlea
import stuff
import traveling_waves as tw

def psth_low_freq():
    bf = 500.0
    fs = 100000.0

    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * bf)
    s = stuff.set_dB_SPL(-8, s)

    ear = cochlea.Sumner2002(hsr=10000, msr=0, lsr=0, freq=500, animal='human')

    hsr, msr, lsr = ear.run(fs, s)

    ax = plt.gca()
    th.plot_psth(hsr['spikes'], bin_size=0.1, axis=ax,
                 histtype='stepfilled',
                 label="PSTH at BF 500 Hz below threshold (-8 dB SPL)")
    ax.legend()
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Spikes per second")


    plt.savefig('psth_at_threshold.eps')
    plt.show()



if __name__ == "__main__":
    psth_low_freq()
