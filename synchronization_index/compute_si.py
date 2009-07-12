import numpy as np
import matplotlib.pyplot as plt

import cochlea
import thorns as th
import stuff
import traveling_waves as tw


def sumner2002_SI(freq_range):
    fs = 100000.0               # Hz
    t = np.arange(0, 0.5, 1.0/fs)

    ear = cochlea.Sumner2002(hsr=1000, msr=0, lsr=0, animal='human')

    dBSPL_range = np.arange(10, 100, 5)

    scores = np.zeros( (len(freq_range), len(dBSPL_range)) )

    for i,freq in enumerate(freq_range):
        for j,dBSPL in enumerate(dBSPL_range):
            print freq, dBSPL,

            # Calculate stimulus
            stim = np.sin(2 * np.pi * t * freq)
            stim = stuff.set_dB_SPL(dBSPL, stim)
            # plt.plot(stim)
            # plt.show()
            #stim = stim * 3000

            ear.bm.set_par("SINGLE_CF", freq)

            # Run ear
            hsr, msr, lsr = ear.run(fs, stim)

            # Compute SI
            si = th.synchronization_index(freq, hsr['spikes'])
            print si

            # Append result
            scores[i,j] = si

    return scores



def holmberg2008_SI():
    fs = 48000.0               # Hz
    t = np.arange(0, 0.5, 1.0/fs) # s

    ear = cochlea.Holmberg2008(hsr=1000, msr=0, lsr=0, animal='human')

    freq_range = tw.real_freq_map
    dBSPL_range = np.arange(10, 100, 5)

    scores = np.zeros( (len(freq_range), len(dBSPL_range)) )

    for i,freq in enumerate(freq_range):
        for j,dBSPL in enumerate(dBSPL_range):
            print freq, dBSPL,

            # Calculate stimulus
            stim = np.sin(2 * np.pi * t * freq)
            stim = stuff.set_dB_SPL(dBSPL, stim)
            # plt.plot(stim)
            # plt.show()
            #stim = stim * 3000

            # Run ear
            hsr, msr, lsr = ear.run(fs, stim, freq_idx=i)

            # Compute SI
            si = th.synchronization_index(freq, hsr['spikes'])
            print si

            # Append result
            scores[i,j] = si

    return scores




if __name__ == "__main__":
    # import cProfile

    freq_range = tw.real_freq_map
    # freq_range = [tw.real_freq_map[61]]

    si_sumner = sumner2002_SI(freq_range)
    si_holmberg = holmberg2008_SI()

    np.savez("si.npz", freq_range=freq_range, si_sumner=si_sumner, si_holmberg=si_holmberg)


    plt.pcolor(si_sumner)
    plt.show()
    plt.pcolor(si_holmberg)
    plt.show()

    scores_max_sumner = np.max(si_sumner, axis=1)
    plt.semilogx(freq_range, scores_max_sumner)
    scores_max_holmberg = np.max(si_holmberg, axis=1)
    plt.semilogx(freq_range, scores_max_holmberg)
    plt.show()

