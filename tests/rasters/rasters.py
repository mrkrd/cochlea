import numpy as np
import matplotlib.pyplot as plt

import cochlea
import thorns as th

def rasters():
    bf = 8000

    fs = 100000.0               # Hz
    t = np.arange(0, 0.025, 1.0/fs)
    tpad = np.arange(0, 0.1, 1.0/fs)
    stim = np.sin(2 * np.pi * t * bf)

    stim = cochlea.set_dB_SPL(50, stim)

    stimpad = np.zeros(tpad.shape)
    stim = np.r_[stimpad, stim, stimpad]

    ear = cochlea.Sumner2002(hsr=1, msr=0, lsr=0, freq=bf,
                             animal='human',
                             sg_type='carney')

    hsr, msr, lsr = ear.run(fs, stim, times=250)

    spikes = hsr['spikes']

    th.plot_raster(spikes)
    th.plot_psth(spikes, bin_size=2)


if __name__ == "__main__":
    rasters()
