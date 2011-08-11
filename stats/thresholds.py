#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import scipy.optimize as opt
import multiprocessing

import thorns as th
import thorns.waves as wv


def calc_threshold_rate(model, model_pars):
    if model.name == "Holmberg2007":
        fs = 48e3
    else:
        fs = 100e3

    ear = model((1000, 0, 0),
                cf=1000,
                **model_pars)

    tmax = 250

    s = np.zeros(fs*tmax/1000)

    anf = ear.run(fs, s)

    rates = [th.calc_rate([train], stimulus_duration=tmax) for train in anf]
    rates = np.array(rates)

    return rates.mean() + rates.std()


def error_function(dbspl, model, model_pars, cf, threshold_rate):
    if model.name == "Holmberg2007":
        fs = 48e3               # Hz
    else:
        fs = 100e3              # Hz

    tone_duration = 250     # ms
    onset = 15                  # ms

    ear = model((1000, 0, 0),
                cf=cf,
                **model_pars)

    s = wv.generate_ramped_tone(fs, cf,
                                tone_duration=tone_duration,
                                ramp_duration=2.5,
                                pad_duration=0,
                                dbspl=dbspl)

    anf = ear.run(fs, s)

    trains = th.trim(anf, onset)
    rate = th.calc_rate(anf, stimulus_duration=(tone_duration-onset))

    error = threshold_rate - rate

    print dbspl, error
    return np.abs(error)


def calc_threshold( (model, freq, model_pars) ):

    threshold_rate = calc_threshold_rate(model, model_pars)

    dbspl_opt = opt.fminbound(error_function,
                              x1=-10,
                              x2=100,
                              args=(model, model_pars, freq, threshold_rate),
                              xtol=0.01)

    return freq, dbspl_opt






def calc_thresholds(model,
                    freqs=np.logspace(np.log10(100), np.log10(10000), 32),
                    model_pars={}):

    pool = multiprocessing.Pool()

    space = [(model, freq, model_pars)
             for freq in freqs]

    thresholds = pool.map(calc_threshold, space)

    thresholds = np.rec.array(thresholds, names=('cfs', 'thresholds'))

    return thresholds


def main():
    import pycat

    model = pycat.Zilany2009
    model_pars = {'powerlaw_implnt':'approx',
                  'with_ffGn':False}

    # threshold_rate = calc_threshold_rate(model, model_pars)

    # print "=== error_function() ==="
    # error_function(dbspl=40,
    #                model=model,
    #                model_pars=model_pars,
    #                cf=1000,
    #                threshold_rate=threshold_rate)

    # print "=== calc_threshold() ==="
    # dbspl_opt = calc_threshold( (model,
    #                              1000,
    #                              model_pars)
    #                         )
    # print "dbspl_opt:", dbspl_opt

    thresholds = calc_thresholds(model=model,
                                 model_pars=model_pars,
                                 freqs=[1000, 3000]
                             )
    print thresholds

if __name__ == "__main__":
    main()
