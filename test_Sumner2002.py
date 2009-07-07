import numpy as np
import matplotlib.pyplot as plt

import ears
import thorns as th

def main():
    fs = 100000.0               # Hz
    t = np.arange(0, 0.1, 1.0/fs)

    freq = 1000

    ear = ears.Sumner2002(hsr=1, msr=0, lsr=0, freq=freq)

    # Calculate stimulus
    stim = np.sin(2 * np.pi * t * freq)
    stim = th.set_dB_SPL(60, stim)

    ear.bm.set_par("SINGLE_CF", freq)

    # Run ear
    hsr, msr, lsr = ear.run(fs, stim, times=2)

    plt.plot(hsr)
    plt.show()

    # Cut PSTH
    psth = hsr[int(len(hsr)*0.3):int(len(hsr)*0.9)]

    # Compute SI
    si = th.calc_SI(fs, freq, psth)
    print si




if __name__ == "__main__":
    main()
