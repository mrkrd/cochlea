# Author:  Marek Rudnicki
# Time-stamp: <2009-11-25 21:25:37 marek>

# Description: Compute the real frequency map of the BM model.

import numpy as np
import matplotlib.pyplot as plt

import traveling_waves as tw
import thorns.waves as w


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
    s = w.set_dB_SPL(80, s)

    s = tw.run_stapes(s)
    s = tw.run_mid_ear_filter(fs, s)
    xBM = tw.run_bm(fs, s, mode='x')

    # plt.imshow(xBM, aspect='auto')
    # plt.show()

    # plt.plot(np.abs(xBM[sec]))
    # plt.show()

    avg = np.max(np.abs(xBM[:,sec]))

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
    plt.plot(real_freq_map)
    plt.show()


if __name__ == "__main__":
    # find_real_freq_map()

    optimize_real_freq_map()
