#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

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


def error_function(


def calc_threshold( (model, freq, model_pars) ):

    threshold_rate = calc_threshold_rate(model, model_pars)

    print threshold_rate






def calc_thresholds(model,
                    freqs=np.logspace(np.log10(100), np.log10(10000), 32),
                    model_pars={}):

    pass




def main():
    import pycat

    _calc_threshold( (pycat.Zilany2009,
                      1000,
                      {'powerlaw_implnt':'approx',
                       'with_ffGn':False}) )

if __name__ == "__main__":
    main()
