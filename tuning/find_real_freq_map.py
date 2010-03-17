# Author:  Marek Rudnicki
# Time-stamp: <2010-03-17 01:43:52 marek>

# Description: Compute the real frequency map of the BM model.

import numpy as np
import biggles

import traveling_waves as tw
import thorns.waves as wv


def calc_mean_displacement(freq, sec, dbspl):
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
    s = tw.set_dbspl(dbspl, s)

    s = tw.scaling_factor * s
    s = tw.run_middle_ear_filter(fs, s)
    xBM = tw.run_bm(fs, s, mode='x')

    avg = np.sum(np.abs(xBM[:,sec]))

    return -avg


def optimize_real_freq_map():
    """
    Runs optimization for each section and finds its BF.
    """
    import scipy.optimize as opt

    dbspl = 0
    real_freq_map = []
    secs = []

    for sec in range(100):
        freq = opt.fminbound(calc_mean_displacement,
                             x1=10,
                             x2=22000,
                             args=(sec,dbspl),
                             xtol=0.1)
        real_freq_map.extend(freq)
        secs.append(sec)

    fname = wv.meta_stamp('maps/real_freq_map.npy', dbspl=dbspl)
    np.save(fname, np.array(real_freq_map))

    p = biggles.FramedPlot()
    p.add( biggles.Curve(secs, real_freq_map) )
    p.show()


if __name__ == "__main__":
    # find_real_freq_map()

    optimize_real_freq_map()
