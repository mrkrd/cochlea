#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import os
import multiprocessing

import thorns.waves as wv
import thorns as th

def _run_model( (model, cf, fs, freq, dbspl, kwargs) ):
    print os.getpid(), freq, dbspl

    tmax = 100
    onset = 10

    ear = model((1500, 0, 0), cf=cf, **kwargs)
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

    return freq, dbspl, rate_hsr


def calc_isointensity_curves(model,
                             cf=3000,
                             fs=100e3,
                             freqs=np.logspace(np.log10(1000), np.log10(8000), 32),
                             dbspls=np.arange(0, 110, 10),
                             **kwargs):

    space = [(model, cf, fs, freq, dbspl, kwargs)
             for m in [model]
             for cf in [cf]
             for fs in [fs]
             for freq in freqs
             for dbspl in dbspls
             for kwargs in [kwargs]]

    pool = multiprocessing.Pool()
    rates = pool.map(_run_model, space)

    rates = np.rec.array(rates, names='freq,dbspl,rate')

    return rates


def main():
    import pycat
    model = pycat.Zilany2009
    pars = {'fs':100e3,
            'powerlaw_implnt':'approx',
            'with_ffGn':False}

    # import cochlea
    # model = cochlea.Sumner2003

    import cochlea
    model = cochlea.Holmberg2007
    pars = {'fs':48000}


    import traveling_waves as tw
    cf = tw.real_freq_map[75]

    # print _run_model( (model, cf, 48000, 3000, 50, {}) )
    # exit()

    rates = calc_isointensity_curves(model, cf, **pars)
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
