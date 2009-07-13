# Author:  Marek Rudnicki
# Time-stamp: <2009-07-06 12:52:35 marek>

# Description: Compute the real frequency map of the BM model.

import numpy as np
import matplotlib.pyplot as plt

import traveling_waves as tw
import stuff

def find_real_freq_map():
    """
    Finds real freq map of BM by scanning range of freq and fining
    which section responds with max displacement.
    """
    fs = 48000.0
    t = np.arange(0, 0.1, 1.0/fs)

    # freq_range = range(1, 22000, 5)
    freq_range = 10**(np.arange(0.1, 4.5, 1e-4))
    # freq_range = [40, 15000]

    scores = np.zeros( (len(freq_range), 100) )

    for freq_idx,freq in enumerate(freq_range):
        print freq

        s = np.sin(2 * np.pi * t * freq)
        s = s * np.hanning(len(s))
        s = stuff.set_dB_SPL(0, s)

        s = tw.run_stapes(s)

        xBM = tw.run_bm(fs, s, mode='x')

        # plt.imshow(xBM, aspect='auto')
        # plt.show()

        avg = np.mean(np.abs(xBM), axis=0)

        # plt.plot(avg)
        # plt.show()

        scores[freq_idx] = avg


    real_freq_map = freq_range[np.argmax(scores, axis=0)]
    np.savez('real_freq_map.npz', scores=scores, freq_range=freq_range, real_freq_map=real_freq_map)



def calc_mean_displacement(freq, sec):
    """
    Computes mean displacement of a given section stimulated by tone
    wiht frequency=freq.

    Note: result is inverted (-) in order to have minium at BF.
    """
    print sec, ":", freq

    fs = 48000.0
    t = np.arange(0, 0.1, 1.0/fs)

    s = np.sin(2 * np.pi * t * freq)
    s = s * np.hanning(len(s))
    s = stuff.set_dB_SPL(0, s)

    s = tw.run_stapes(s)

    xBM = tw.run_bm(fs, s, mode='x')

    # plt.imshow(xBM, aspect='auto')
    # plt.show()

    # plt.plot(np.abs(xBM[sec]))
    # plt.show()

    avg = np.mean(np.abs(xBM[:,sec]))

    return -avg


def optimize_real_freq_map():
    """
    Runs optimization for each section and finds its BF.
    """
    import scipy.optimize as opt

    sec_num = 100
    last_x0 = 1000
    real_freq_map = np.zeros(sec_num)
    for sec in range(sec_num):
        real_freq_map[sec] = opt.fminbound(calc_mean_displacement,
                                           x1=10,
                                           x2=22000,
                                           args=(sec,),
                                           xtol=0.01)

    np.save('real_freq_map.npy', real_freq_map)


if __name__ == "__main__":
    # find_real_freq_map()

    optimize_real_freq_map()
