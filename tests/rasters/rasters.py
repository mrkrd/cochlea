import numpy as np

import cochlea
import thorns as th
import stuff

def rasters():
    bf = 500

    fs = 100000.0               # Hz
    t = np.arange(0, 0.3, 1.0/fs)
    tpad = np.arange(0, 0.1, 1.0/fs)
    stim = np.sin(2 * np.pi * t * bf)

    stim = stuff.set_dB_SPL(50, stim)

    stimpad = np.zeros(tpad.shape)
    stim = np.r_[stimpad, stim, stimpad]

    ear = cochlea.Sumner2002(hsr=1, msr=0, lsr=0, freq=bf)

    hsr, msr, lsr = ear.run(fs, stim, times=250)

    spikes = hsr['spikes']

    th.plot_raster(spikes)
    th.plot_psth(spikes, bin_size=2)


if __name__ == "__main__":
    rasters()
