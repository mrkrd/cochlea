import numpy as np

import gp_ear
import thorns as th


def rasters():
    bf = 1000
    trial_num = 1

    fs = 100000.0               # Hz
    t = np.arange(0, 0.3, 1.0/fs)
    tpad = np.arange(0, 0.1, 1.0/fs)
    stim = np.sin(2 * np.pi * t * bf)
    stim = th.set_dB_SPL(0, stim)
    stimpad = np.zeros(tpad.shape)
    stim = np.r_[stimpad, stim, stimpad]

    ear = gp_ear.Sumner2002(hsr=250, msr=0, lsr=0, freq=bf)

    hsr_list = []

    for i in range(trial_num):
        hsr, msr, lsr = ear.run(fs, stim)

        hsr_list.append(hsr)


    spikes = th.timef_to_spikes(fs, np.asarray(hsr_list).astype(int))
    # th.plot_raster(spikes)
    th.plot_psth(spikes, bin_size=2)


if __name__ == "__main__":
    rasters()
