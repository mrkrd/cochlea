import numpy as np
import matplotlib.pyplot as plt

import cochlea
import thorns as th

def main():
    fs = 100000.0               # Hz
    t = np.arange(0, 0.1, 1.0/fs)

    freq = 1000

    ear = cochlea.Sumner2002(hsr=1, msr=0, lsr=0, freq=freq)

    # Calculate stimulus
    stim = np.sin(2 * np.pi * t * freq)
    stim = cochlea.set_dB_SPL(60, stim)

    # Run ear
    hsr, msr, lsr = ear.run(fs, stim, times=250)


    th.plot_raster(hsr['spikes'])
    th.plot_psth(hsr['spikes'])



if __name__ == "__main__":
    main()
