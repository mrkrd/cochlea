#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import scipy.optimize as opt
import multiprocessing

import thorns as th
import thorns.waves as wv
from . import binary

def calc_threshold_rate(model, model_pars):
    if model.name == "Holmberg2007":
        import traveling_waves as tw
        fs = 48e3
        cf = tw.find_closest_freq_in_map(1000)
    else:
        fs = 100e3
        cf = 1000

    ear = model((10000, 0, 0),
                cf=cf,
                **model_pars)

    tmax = 250

    s = np.zeros(fs*tmax/1000)

    anf = ear.run(s, fs)

    trains = anf['spikes']
    durations = anf['duration']

    rates = [1000*len(train)/duration
             for train,duration in zip(trains,durations)]
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

    anf = ear.run(s, fs)

    trains = th.trim(anf, onset)
    rate = th.stats.rate(anf)

    error = rate - threshold_rate

    print dbspl, error
    return error


def calc_threshold( (model, freq, model_pars) ):

    threshold_rate = calc_threshold_rate(model, model_pars)
    print "threshold_rate:", threshold_rate

    dbspl_opt = binary.find_threshold(error_function,
                                      args=(model, model_pars, freq, threshold_rate),
                                      init_range=(-100, 100),
                                      desired_range=0.1)


    return freq, dbspl_opt






def calc_thresholds(model,
                    freqs=np.logspace(np.log10(100), np.log10(16000), 32),
                    model_pars={}):

    pool = multiprocessing.Pool()

    space = [(model, freq, model_pars)
             for freq in freqs]

    thresholds = pool.map(calc_threshold, space)

    thresholds = np.rec.array(thresholds, names=('cfs', 'thresholds'))

    return thresholds



def calc_human_hearing_threshold(freqs):
    """Terhardt, E. (1979). Calculating virtual pitch. Hearing
    Research, 1(2):155-182.

    http://www.diracdelta.co.uk/science/source/t/h/threshold%20of%20hearing/source.html

    """
    f = freqs / 1000                # kHz -> Hz

    th = 3.64 * f**(-0.8) - \
         6.5 * np.exp(-0.6 * (f - 3.3)**2) + \
         10**(-3) * f**4

    return th



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

    print "=== calc_threshold() ==="
    dbspl_opt = calc_threshold( (model,
                                 10000,
                                 model_pars)
                            )
    print "dbspl_opt:", dbspl_opt

    # thresholds = calc_thresholds(model=model,
    #                              model_pars=model_pars,
    #                              freqs=[1000, 3000]
    #                          )
    # print thresholds

if __name__ == "__main__":
    main()
