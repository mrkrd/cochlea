#!/usr/bin/env python

from __future__ import division
from __future__ import print_function


raise NotImplementedError


__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd
import logging

import thorns as th
import thorns.waves as wv



def calc_spont_threshold(model, model_pars=None):

    if model_pars is None:
        model_pars = {}

    pars = dict(model_pars)

    fs = pars.setdefault('fs', 100e3)
    cf = pars.setdefault('cf', 1e3)

    tmax = 250e-3
    s = np.zeros(fs*tmax)

    anf = model(
        sound=s,
        anf_num=(10000,0,0),
        seed=0,
        **pars
    )

    sis = [th.calc_si(train,cf) for _,train in anf.iterrows()]
    sis = np.array(sis)

    threshold = sis.mean() + sis.std()

    return threshold



def error_func(dbspl, model, cf, model_pars, spont_si):

    if model_pars is None:
        model_pars = {}

    pars = dict(model_pars)

    fs = pars.setdefault('fs', 100e3)

    tone_duration = 250e-3
    onset = 15e-3

    s = wv.make_ramped_tone(
        fs,
        cf,
        duration=tone_duration,
        ramp=2.5e-3,
        pad=0,
        dbspl=dbspl
    )

    anf = model(
        sound=s,
        anf_num=(1000, 0, 0),
        cf=cf,
        seed=0,
        **pars
    )

    trains = th.trim(anf, onset, None)
    si = th.calc_si(trains, cf)

    error = si - spont_si

    logging.debug("{} {} {}".format(dbspl, si, error))

    return error



def calc_hearing_threshold_si(model, cf, spont_si, model_pars=None):

    kwargs = {
        'model':model,
        'model_pars': model_pars,
        'cf': cf,
        'spont_si': spont_si
    }
    dbspl_opt = th.util.find_zero(
        error_func,
        kwargs=kwargs,
        x1=-100,
        x2=50,
        xtol=0.1
    )

    return dbspl_opt






def calc_hearing_thresholds_si(
        model,
        cfs=None,
        **model_pars
):

    if cfs is None:
        cfs = np.logspace(np.log10(100), np.log10(16000), 32)


    spont_si = th.util.apply(
        calc_spont_threshold,
        model=model,
        model_pars=model_pars
    )

    space = [
        {
            'model': model,
            'model_pars': model_pars,
            'spont_si': spont_si,
            'cf': cf
        }
        for cf in cfs
    ]


    thresholds = th.util.map(
        calc_hearing_threshold_si,
        space,
        # backend='multiprocessing'
    )
    thresholds = list(thresholds)

    thresholds = pd.Series(thresholds, index=cfs)

    return thresholds





def main():
    import cochlea
    import matplotlib.pyplot as plt

    spont_si = calc_spont_threshold(
        cochlea.run_zilany2009
    )
    print(spont_si)


    r = calc_hearing_threshold_si(
        model=cochlea.run_zilany2009,
        cf=4e3,
        spont_si=spont_si
    )

    print(r)
    exit()


    ths = calc_hearing_thresholds_si(
        model=cochlea.run_zilany2009
    )
    print(ths)

    ths.plot(logx=True)
    plt.show()



if __name__ == "__main__":
    main()
