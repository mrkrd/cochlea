#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd

import marlib as mr
import marlib.thorns as th
import marlib.waves as wv

def _run_model(model, dbspl, cf, fs):

    duration = 100e-3
    onset = 10e-3

    s = wv.make_ramped_tone(
        fs=fs,
        freq=cf,
        duration=duration,
        pad=0,
        dbspl=dbspl
    )

    anf = model(
        sound=s,
        fs=fs,
        anf_num=(250,250,250),
        cf=cf,
        seed=0
    )

    hsr = anf[anf['type']=='hsr']
    hsr = th.trim(hsr, onset, None)
    rate_hsr = th.calc_rate(hsr)

    msr = anf[anf['type']=='msr']
    msr = th.trim(msr, onset, None)
    rate_msr = th.calc_rate(msr)

    lsr = anf[anf['type']=='lsr']
    lsr = th.trim(lsr, onset, None)
    rate_lsr = th.calc_rate(lsr)

    out = {
        'cf': cf,
        'dbspl': dbspl,
        'hsr': rate_hsr,
        'msr': rate_msr,
        'lsr': rate_lsr
    }
    return out


def calc_rate_intensity(model, cfs, fs=100e3, dbspls=None):

    if dbspls is None:
        dbspls = np.arange(-10, 100, 5)

    space = [
        {
            'model': model,
            'dbspl': dbspl,
            'cf': cf,
            'fs': fs,
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

    rates = calc_rate_intensity(
        model=cochlea.run_holmberg2007,
        cfs=[cochlea.holmberg2007.real_freq_map[80]],
        fs=48e3,
        # dbspls=[0, 20, 50]
    )


    rates.hsr.plot()
    rates.msr.plot()
    rates.lsr.plot()
    mr.show()

    print(rates)

if __name__ == "__main__":
    main()
