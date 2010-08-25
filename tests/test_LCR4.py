from __future__ import division

import numpy as np
import biggles

from traveling_waves import _tw, S_ST, bm_pars

def main():
    fs = 48000                  # Hz

    forward = np.loadtxt('forward.txt') * S_ST
    t = np.arange(len(forward)) / fs * 1000

    ref = np.loadtxt('LCR4_70.txt')


    xbm = _tw.bm_wave(fs=48000,
                      signal=forward)

    xbm = xbm[:,70]

    lcr4 = _tw._single_LCR4(fs, xbm, 70)

    sel = range(400, 2000)
    t = t[sel]
    lcr4 = lcr4[sel]
    ref = ref[sel]
    p = biggles.FramedPlot()
    p.add( biggles.Curve(t, ref, color='red', width=3) )
    p.add( biggles.Curve(t, lcr4))
    p.show()




if __name__ == "__main__":
    main()
