#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import scipy as sp
import scipy.signal
import matplotlib.pyplot as plt

import traveling_waves as tw

def calc_curve(fs):
    b, a = tw._calc_middle_ear_coefs(fs)

    h, w = sp.signal.freqz(b, a)

    freq = fs * h / (2*np.pi)
    amp = np.abs(w)

    plt.plot(freq, amp)



def main():
    calc_curve(48000)
    calc_curve(100000)

    plt.show()

if __name__ == "__main__":
    main()
