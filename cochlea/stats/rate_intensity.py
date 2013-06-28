#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd

import mrlib as mr
import mrlib.thorns as th
import mrlib.waves as wv

def _run_model(model, dbspl, cf, model_pars):

    duration = 100e-3
    onset = 10e-3

    fs = model_pars.setdefault('fs', 100e3)
    model_pars.setdefault('anf_num', (250,250,250))
    model_pars.setdefault('seed', 0)

    s = wv.ramped_tone(
        fs=fs,
        freq=cf,
        duration=duration,
        pad=0,
        dbspl=dbspl
    )

    anf = model(
        sound=s,
        cf=cf,
        **model_pars
    )

    hsr = anf[anf['type']=='hsr']
    hsr = th.trim(hsr, onset, None)
    rate_hsr = th.rate(hsr)

    msr = anf[anf['type']=='msr']
    msr = th.trim(msr, onset, None)
    rate_msr = th.rate(msr)

    lsr = anf[anf['type']=='lsr']
    lsr = th.trim(lsr, onset, None)
    rate_lsr = th.rate(lsr)

    out = {
        'cf': cf,
        'dbspl': dbspl,
        'hsr': rate_hsr,
        'msr': rate_msr,
        'lsr': rate_lsr
    }
    return out


def calc_rate_intensity(model, cfs, dbspls=None, **model_pars):

    if dbspls is None:
        dbspls = np.arange(-10, 100, 5)

    space = [
        {
            'model': model,
            'dbspl': dbspl,
            'cf': cf,
            'model_pars': model_pars,
        }
        for dbspl in dbspls
        for cf in cfs
    ]

    rates = mr.map(
        _run_model,
        space
    )

    rates = pd.DataFrame(list(rates))
    return rates


def main():
    import cochlea

    # rates = calc_rate_intensity(
    #     model=cochlea.run_holmberg2007,
    #     cfs=[cochlea.holmberg2007.real_freq_map[80]],
    #     fs=48e3,
    #     # dbspls=[0, 20, 50]
    # )

    rates = calc_rate_intensity(
        model=cochlea.run_zilany2013,
        cfs=[5e3],
        fs=100e3,
        species='human'
        # dbspls=[0, 20, 50]
    )


    rates.hsr.plot()
    rates.msr.plot()
    rates.lsr.plot()
    mr.show()

    print(rates)

if __name__ == "__main__":
    main()
