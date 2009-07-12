import numpy as np
import matplotlib.pyplot as plt

import cochlea
import thorns as th
import stuff
import traveling_waves as tw

def thresholds_Sumner2002(freq_range):
    fs = 100000.0

    hsr_num = 1000
    ear = cochlea.Sumner2002(hsr=hsr_num, msr=0, lsr=0, freq=1000, animal='human')

    # Measure spontanious activity
    tmax = 0.5
    t = np.arange(0, tmax, 1/fs)
    s = np.zeros_like(t)
    hsr, msr, lsr = ear.run(fs, s)
    all_spikes = np.concatenate(tuple(hsr['spikes']))
    hist = np.histogram(all_spikes, bins=np.ceil(tmax * 1000))[0]
    mean = hist.mean()
    sd = hist.std()
    sps_threshold = mean + sd

    print mean, sd, sps_threshold, len(all_spikes)


    # Find dB_SPL for which threshold is exceded
    # freq_range = [1000]
    # freq_range = range(100, 20000, 200)

    dB_SPL_threshold = 25
    dB_SPL_step = 0.5

    thresholds = []
    for freq in freq_range:
        ear.set_freq(freq)
        s = np.sin(2 * np.pi * t * freq)

        dB_SPL = dB_SPL_threshold
        no_change = True
        trend = None

        print freq, ":",

        while no_change:

            # Compute spikes per second at current dB_SPL
            s = stuff.set_dB_SPL(dB_SPL, s)
            hsr, msr, lsr = ear.run(fs, s)
            spike_num = np.concatenate(tuple(hsr['spikes'])).size
            sps = spike_num / (tmax * 1000)

            print dB_SPL,

            # Check the trend (up/down)
            if sps > sps_threshold:
                if trend == 'up':
                    no_change = False
                    dB_SPL_threshold = dB_SPL
                dB_SPL -= dB_SPL_step
                trend = 'down'
            else:
                dB_SPL += dB_SPL_step
                if trend == 'down':
                    no_change = False
                    dB_SPL_threshold = dB_SPL
                trend = 'up'


        print "=", dB_SPL_threshold
        thresholds.append(dB_SPL_threshold)

    return thresholds




def thresholds_Holmberg2008(freq_range):
    fs = 48000.0

    hsr_num = 1000
    ear = cochlea.Holmberg2008(hsr=hsr_num, msr=0, lsr=0, animal='human')

    # Measure spontanious activity
    tmax = 0.5
    t = np.arange(0, tmax, 1/fs)
    s = np.zeros_like(t)
    ear.set_freq(tw.bm_pars.real_freq_map[10])
    hsr, msr, lsr = ear.run(fs, s)
    all_spikes = np.concatenate(tuple( hsr['spikes'] ))
    hist = np.histogram(all_spikes, bins=np.ceil(tmax*1000))[0]
    mean = hist.mean()
    sd = hist.std()
    sps_threshold = mean + sd

    print mean, sd, sps_threshold, len(all_spikes)


    # Find dB_SPL for which threshold is exceded
    #freq_range = [1000]
    # freq_range = range(100, 20000, 200)

    dB_SPL_threshold = -4
    dB_SPL_step = 0.5

    thresholds = []
    for freq in freq_range:
        s = np.sin(2 * np.pi * t * freq)
        ear.set_freq(freq)

        dB_SPL = dB_SPL_threshold
        no_change = True
        trend = None

        print freq, ":",

        while no_change:

            # Compute spikes per second at current dB_SPL
            s = stuff.set_dB_SPL(dB_SPL, s)
            hsr, msr, lsr = ear.run(fs, s)

            spike_num = np.concatenate(tuple(hsr['spikes'])).size
            sps = spike_num / (tmax*1000)

            print dB_SPL

            # Check the trend (up/down)
            if sps > sps_threshold:
                if trend == 'up':
                    no_change = False
                    dB_SPL_threshold = dB_SPL
                dB_SPL -= dB_SPL_step
                trend = 'down'
            else:
                dB_SPL += dB_SPL_step
                if trend == 'down':
                    no_change = False
                    dB_SPL_threshold = dB_SPL
                trend = 'up'


        print "=", dB_SPL_threshold
        thresholds.append(dB_SPL_threshold)

    return thresholds




if __name__ == "__main__":
    freq_range = tw.bm_pars.real_freq_map
    # Note: index 38 -> ~1000Hz [990Hz]
    # freq_range = [tw.bm_pars.real_freq_map[38]]
    print freq_range

    sumner_th = thresholds_Sumner2002(freq_range)
    holmberg_th = thresholds_Holmberg2008(freq_range)

    np.savez('th.npz', sumner_th=sumner_th, holmberg_th=holmberg_th)

    plt.semilogx(freq_range, sumner_th)
    plt.semilogx(freq_range, holmberg_th)
    plt.show()

