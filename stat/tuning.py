#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np


def _calc_point( (model, cf, freq, mean_rate, sd_rate, kwargs) ):
    print os.getpid()

    tmax = 250                  # ms
    onset = 30                  # ms
    step = 1                    # dB
    threshold_rate = mean_rate + sd_rate
    trend = None
    no_change = True

    ear = model((1000, 0, 0), cf=cf, **kwargs)


    while no_change:
        # Compute spikes per second at current dB SPL
        s = wv.generate_ramped_tone(fs, freq,
                                    tone_duration=tmax,
                                    pad_duration=0,
                                    dbspl=dbspl)
        anf = ear.run(fs, s)
        anf = th.trim(anf, onset)
        rate = th.calc_rate(anf, stimulus_duration=(tmax-onset))


        # Check the trend (up/down)
        if rate > threshold_rate:
            if trend == 'up':
                no_change = False
                threshold = dbspl
            dbspl -= step
            trend = 'down'
        else:
            dbspl += step
            if trend == 'down':
                no_change = False
                threshold = dbspl
            trend = 'up'


    return freq, threshold


def _calc_spont_rate(model, **kwargs):


    ear = model((1000, 0, 0), **kwargs)


def main():
    import pycat

    ear =


if __name__ == "__main__":
    main()
