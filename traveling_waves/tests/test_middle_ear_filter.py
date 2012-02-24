#!/usr/bin/env python

from __future__ import division
__author__ = "Marek Rudnicki"

import numpy as np
import biggles

from traveling_waves import run_middle_ear_filter, _tw, S_ED

def main():
    fs = 48000

    signal = np.loadtxt('signal.txt')
    t = np.arange(len(signal)) / fs * 1000

    ref = np.loadtxt('forward.txt')

    forward = run_middle_ear_filter(fs, signal) * S_ED

    sec = range(450,530)
    t = t[sec]
    forward = forward[sec]
    ref = ref[sec]
    p = biggles.FramedPlot()
    p.add( biggles.Curve(t, ref, color='red', width=3))
    p.add( biggles.Curve(t, forward))
    p.show()

if __name__ == "__main__":
    main()
