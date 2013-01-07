#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd

import marlib.thorns as th
import marlib.waves as wv

from binary import find_zero

def calc_spont_threshold(model, model_pars):
    if 'fs' not in model_pars:
        model_pars['fs'] = 100e3

    if 'cf' not in model_pars:
        model_pars['cf'] = 1e3

    if 'seed' not in model_pars:
        model_pars['seed'] = 0


    tmax = 250e-3
    s = np.zeros(model_pars['fs']*tmax)

    anf = model(
        sound=s,
        anf_num=(10000, 0, 0),
        **model_pars
    )

    rates = [th.calc_rate(train) for _,train in anf.iterrows()]
    rates = np.array(rates)

    threshold = rates.mean() + rates.std()

    return threshold



def error_func(dbspl, model, cf, model_pars, spont_rate):

    if 'fs' not in model_pars:
        model_pars['fs'] = 100e3

    if 'seed' not in model_pars:
        model_pars['seed'] = 0

    fs = model_pars['fs']

    tone_duration = 250e-3
    onset = 15e-3

    s = wv.make_ramped_tone(
        fs, cf,
        duration=tone_duration,
        ramp=2.5e-3,
        pad=0,
        dbspl=dbspl
    )

    anf = model(
        sound=s,
        anf_num=(1000, 0, 0),
        cf=cf,
        **model_pars
    )

    trains = th.trim(anf, onset, None)
    rate = th.calc_rate(trains)

    error = rate - spont_rate

    print(dbspl, rate, error)

    return error



def calc_hearing_threshold_rate(model, cf, model_pars, spont_rate):

    kwargs = {
        'model':model,
        'model_pars': model_pars,
        'cf': cf,
        'spont_rate': spont_rate
    }
    dbspl_opt = find_zero(
        error_func,
        kwargs=kwargs,
        x1=-100,
        x2=100,
        xtol=0.1
    )

    return dbspl_opt






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
    import cochlea

    spont_rate = calc_spont_threshold(
        cochlea.run_zilany2009,
        {}
    )
    print(spont_rate)


    r = calc_hearing_threshold_rate(
        model=cochlea.run_zilany2009,
        cf=1e3,
        model_pars={},
        spont_rate=spont_rate
    )

    print(r)

    exit()

    # print "=== error_function() ==="
    # error_function(dbspl=40,
    #                model=model,
    #                model_pars=model_pars,
    #                cf=1000,
    #                threshold_rate=threshold_rate)

    print("=== calc_threshold() ===")
    dbspl_opt = calc_threshold(
        (
            model,
            10000,
            model_pars
        )
    )
    print("dbspl_opt:", dbspl_opt)

    # thresholds = calc_thresholds(model=model,
    #                              model_pars=model_pars,
    #                              freqs=[1000, 3000]
    #                          )
    # print thresholds

if __name__ == "__main__":
    main()
