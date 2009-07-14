# Author:  Marek Rudnicki
# Time-stamp: <2009-07-14 22:53:19 marek>

# Description: Compute center frequencies of LCR module.

import numpy as np
import matplotlib.pyplot as plt

import traveling_waves as tw
import tw._bm
import stuff


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

    xBM = tw.run_bm(fs, s, mode='x', with_LCR=False)

    # plt.imshow(xBM, aspect='auto')
    # plt.show()

    # plt.plot(np.abs(xBM[sec]))
    # plt.show()

    avg = np.mean(np.abs(xBM[:,sec]))

    return -avg


def optimize_freq_map():
    """
    Runs optimization for each section and finds its BF.  LCR
    resonators are *disabled*.
    """
    import scipy.optimize as opt

    sec_num = 100
    last_x0 = 1000
    freq_map = np.zeros(sec_num)
    for sec in range(sec_num):
        freq_map[sec] = opt.fminbound(calc_mean_displacement,
                                      x1=10,
                                      x2=22000,
                                      args=(sec,),
                                      xtol=0.01)

    np.save('freq_map_without_LCR.npy', freq_map)


def plot_freq_map():
    freq_map = np.load('freq_map_without_LCR.npy')

    ax = plt.gca()

    ax.plot(tw.real_freq_map, label="LCR on")
    ax.plot(freq_map, label="LCR off")
    ax.set_ylabel("Frequency [Hz]")
    ax.set_xlabel("Section number")
    ax.legend(loc='best')

    plt.savefig("freq_map_without_LCR.eps")
    plt.show()


if __name__ == "__main__":
    # find_real_freq_map()

    optimize_freq_map()

    plot_freq_map()
