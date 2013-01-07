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

    cf = model_pars['cf']

    tmax = 250e-3
    s = np.zeros(model_pars['fs']*tmax)

    anf = model(
        sound=s,
        anf_num=(10000, 0, 0),
        **model_pars
    )

    sis = [th.calc_si(train,cf) for _,train in anf.iterrows()]
    sis = np.array(sis)

    threshold = sis.mean() + sis.std()

    return threshold



def error_func(dbspl, model, cf, model_pars, spont_si):

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
    si = th.calc_si(trains, cf)

    error = si - spont_si

    print(dbspl, si, error)

    return error



def calc_hearing_threshold_si(model, cf, model_pars, spont_si):

    kwargs = {
        'model':model,
        'model_pars': model_pars,
        'cf': cf,
        'spont_si': spont_si
    }
    dbspl_opt = find_zero(
        error_func,
        kwargs=kwargs,
        x1=-100,
        x2=70,
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






def main():
    import cochlea

    spont_si = calc_spont_threshold(
        cochlea.run_zilany2009,
        {}
    )
    print(spont_si)


    r = calc_hearing_threshold_si(
        model=cochlea.run_zilany2009,
        cf=1e3,
        model_pars={},
        spont_si=spont_si
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
