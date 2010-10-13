#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import os
import multiprocessing

import thorns as th
import thorns.waves as wv


def _calc_point( (model, fs, cf, freq, mean_rate, sd_rate, pars) ):
    print os.getpid()

    tmax = 250                  # ms
    onset = 30                  # ms
    step = 1                    # dB
    threshold_rate = mean_rate + sd_rate
    trend = None
    no_change = True
    dbspl = 0
    ear = model((1000, 0, 0), cf=cf, **pars)


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


def _calc_spont_rate(model, fs, pars):

    ear = model((1000, 0, 0), cf=1000, **pars)

    tmax = 250

    s = np.zeros(fs*tmax/1000)

    anf = ear.run(fs, s)

    rates = [th.calc_rate([train], stimulus_duration=tmax) for train in anf]
    rates = np.array(rates)

    return rates.mean(), rates.std()



def calc_tuning(model,
                cf,
                fs=100e3,
                freqs=np.logspace(np.log10(500), np.log10(6000), 32),
                **pars):

    mean,sd = _calc_spont_rate(model, fs=fs, pars=pars)

    space = [(model, fs, cf, freq, mean_rate, sd_rate, pars)
             for m in [model]
             for fs in [fs]
             for cf in [cf]
             for freq in freqs
             for mean_rate in [mean]
             for sd_rate in [sd]
             for pars in [pars]]

    pool = multiprocessing.Pool()
    thresholds = pool.map(_calc_point, space)

    thresholds = np.rec.array(thresholds, names='freq,threshold')

    return thresholds



def main():
    import pycat
    model = pycat.Zilany2009
    fs = 100e3
    pars = { 'powerlaw_implnt':'approx',
             'with_ffGn':False }
    cf = 4000

    mean,sd = _calc_spont_rate(model, fs=fs, pars=pars)
    print _calc_point( (model, fs, cf, 3500, mean, sd, pars) )


    print '==========================================='
    print calc_tuning(model, cf=cf, fs=fs,
                      freqs=[3500, 4000],
                      **pars)



if __name__ == "__main__":
    main()
