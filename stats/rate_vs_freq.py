#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import os
import multiprocessing

import thorns.waves as wv
import thorns as th

def _run_model( (model, cf, freq, dbspl) ):
    print os.getpid(), freq, dbspl

    tmax = 100
    fs = 48000
    onset = 10

    ear = model((1500, 0, 0), cf=cf)
    # ear = model((250, 0, 0), cf=cf,
    #             powerlaw_implnt='approx',
    #             with_ffGn=False)

    s = wv.generate_ramped_tone(fs, freq,
                                tone_duration=tmax,
                                pad_duration=0,
                                dbspl=dbspl)
    anf = ear.run(fs, s)

    hsr = anf.where(typ='hsr')
    hsr = th.trim(hsr, onset)
    rate_hsr = th.calc_rate(hsr, stimulus_duration=(tmax-onset))

    return rate_hsr


def calc_rate_curves(model, cf):
    space = [(model, cf, freq, dbspl)
             for m in [model]
             for cf in [cf]
             for freq in np.logspace(np.log10(5000), np.log10(6000), 32)
             for dbspl in np.arange(0, 90, 20)]

    pool = multiprocessing.Pool()
    rates = pool.map(_run_model, space)

    freqs = [el[2] for el in space]
    dbspls = [el[3] for el in space]
    rates = zip(freqs,dbspls,rates)
    rates = np.rec.array(rates, names='freq,dbspl,rate')

    return rates


def main():
    # import pycat
    # model = pycat.Zilany2009

    # import cochlea
    # model = cochlea.Sumner2003

    import cochlea
    model = cochlea.Holmberg2007
    import traveling_waves
    cf = traveling_waves.real_freq_map[75]

    # print _run_model( (model, 5000, 3000, 50) )
    # exit()

    rates = calc_rate_curves(model, cf)
    print rates

    import biggles
    p = biggles.FramedPlot()
    p.xlog = 1
    for dbspl in np.unique(rates['dbspl']):
        selected = rates[ rates['dbspl']==dbspl ]
        p.add( biggles.Curve(selected['freq'], selected['rate']) )

    p.show()

if __name__ == "__main__":
    main()
