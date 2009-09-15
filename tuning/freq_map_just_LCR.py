# Author:  Marek Rudnicki
# Time-stamp: <2009-09-15 15:37:01 marek>

# Description: Compute center frequencies of LCR module.

import numpy as np
import matplotlib.pyplot as plt

import traveling_waves as tw
from traveling_waves import bm_pars, _bm
import stuff


def calc_mean_displacement(freq, sec):
    """
    Computes mean displacement of a given section stimulated by tone
    wiht frequency=freq.

    Note: result is negated in order to have minium (instead of max) at BF.
    """
    print sec, ":", freq

    fs = 48000.0
    t = np.arange(0, 0.1, 1.0/fs)

    s = np.sin(2 * np.pi * t * freq)
    s = s * np.hanning(len(s))
    s = stuff.set_dB_SPL(0, s)
    s = tw.run_stapes(s)

    signal = np.zeros( (s.size, 100) )
    signal[:,sec] = s

    _bm.LCR4_init(fs,
                  bm_pars.freq_map,
                  bm_pars.Qmin,
                  bm_pars.SAT1,
                  bm_pars.SAT4)


    signal = _bm.LCR4(signal,
                      bm_pars.Qmax,
                      bm_pars.Qmin)

    # plt.imshow(xBM, aspect='auto')
    # plt.show()

    # plt.plot(signal[:,sec])
    # plt.plot(out[:,sec])
    # plt.show()

    avg = np.mean(np.abs(signal[:,sec]))

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

    np.save('freq_map_just_LCR.npy', freq_map)


def plot_freq_map():
    freq_map = np.load('freq_map_just_LCR.npy')

    ax = plt.gca()

    ax.plot(tw.real_freq_map[::-1], label="real_freq_map: BM with LCR")
    ax.plot(freq_map, label="Only LCR")
    ax.plot(tw.bm_pars.freq_map, label="freq_map: LCR4 input parameter")
    ax.plot(1/tw.bm_pars.delay_time[::-1], label="1 / delay")
    ax.set_ylabel("Frequency [Hz]")
    ax.set_xlabel("Section number")
    ax.legend(loc='best')

    plt.savefig("freq_map_just_LCR.eps")
    plt.show()


def plot_single_sec_characteristic():

    sec = 20

    characteristic = []
    for freq in tw.real_freq_map:
        characteristic.append(calc_mean_displacement(freq, sec))

    ax = plt.gca()
    ax.semilogx(tw.real_freq_map, characteristic,
                label=str(tw.real_freq_map[99-sec]))
    ax.legend(loc='best')

    plt.show()

if __name__ == "__main__":
    # find_real_freq_map()

    # optimize_freq_map()

    plot_freq_map()

    # plot_single_sec_characteristic()
