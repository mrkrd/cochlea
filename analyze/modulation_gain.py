#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import multiprocessing
import os

import thorns as th
import thorns.waves as wv


def _run_model( (model, fm) ):
    print os.getpid(), fm

    tmax = 600                  # ms
    cf = 10000                  # Hz
    fs = 1e5                    # Hz
    m = 1                       # modulation depth

    ear = model((200, 0, 0), cf=cf,
                powerlaw_implnt='approx',
                with_ffGn=False)

    s = wv.generate_am_tone(fs,
                            fm=fm,
                            fc=cf,
                            modulation_depth=m,
                            tone_duration=tmax)
    s = wv.set_dbspl(10, s)
    anf = ear.run(fs, s)

    si = th.calc_si(anf['spikes'], fm)
    gain = 20 * np.log10(2*si / m)

    return gain


def calc_modulation_gain(model):
    fms = np.logspace(1, 3.3, 16)

    space = [(m, fm)
             for m in [model]
             for fm in fms]

    pool = multiprocessing.Pool()

    gains = pool.map(_run_model, space)

    gains = [(fm,gain) for fm,gain in zip(fms,gains)]
    gains = np.rec.array(gains, names='fm,gain')

    return gains


def main():
    import biggles

    # import cochlea
    # gains = calc_modulation_gain(cochlea.Sumner2003)

    import pycat
    gains = calc_modulation_gain(pycat.Zilany2009)

    plot = biggles.FramedPlot()
    plot.add( biggles.Curve(gains['fm'], gains['gain']) )
    plot.xlog = 1
    plot.xlabel = "Modulation Frequency (Hz)"
    plot.ylabel = "Modulation Gain (dB)"

    plot.write_eps("plots/modulation_gain.eps")
    plot.show()

if __name__ == "__main__":
    main()
