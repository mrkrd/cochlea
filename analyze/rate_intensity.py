#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import multiprocessing
import os

import thorns as th
import thorns.waves as wv

def _run_model( (model, dbspl) ):
    print os.getpid(), dbspl

    tmax = 100
    cf = 8000
    fs = 100000
    onset = 10

    ear = model((250, 250, 250), cf=cf)

    s = wv.generate_ramped_tone(fs, cf,
                                tone_duration=tmax,
                                pad_duration=0,
                                dbspl=dbspl)
    anf = ear.run(fs, s)

    hsr = anf[anf['typ']=='hsr']['spikes']
    hsr = th.trim(hsr, onset)
    rate_hsr = th.calc_rate(hsr, stimulus_duration=(tmax-onset))

    msr = anf[anf['typ']=='msr']['spikes']
    msr = th.trim(msr, onset)
    rate_msr = th.calc_rate(msr, stimulus_duration=(tmax-onset))

    lsr = anf[anf['typ']=='lsr']['spikes']
    lsr = th.trim(lsr, onset)
    rate_lsr = th.calc_rate(lsr, stimulus_duration=(tmax-onset))

    return rate_hsr, rate_msr, rate_lsr


def calc_rate_intensity(model):
    dbspls = np.arange(-10, 100, 5)

    space = [(m, dbspl)
             for m in [model]
             for dbspl in dbspls]

    pool = multiprocessing.Pool()
    rates = pool.map(_run_model, space)

    rates = [(dbspl,)+rate for dbspl,rate in zip(dbspls,rates)]
    rates = np.rec.array(rates, names='dbspl,hsr,msr,lsr')

    return rates


def main():
    import biggles
    import cochlea

    rates = calc_rate_intensity(cochlea.Sumner2002)

    plot = biggles.FramedPlot()
    plot.add( biggles.Curve(rates['dbspl'], rates['hsr']) )
    plot.add( biggles.Curve(rates['dbspl'], rates['msr']) )
    plot.add( biggles.Curve(rates['dbspl'], rates['lsr']) )

    plot.show()

if __name__ == "__main__":
    main()
