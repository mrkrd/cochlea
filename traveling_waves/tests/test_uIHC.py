from __future__ import division

import numpy as np
import biggles

from traveling_waves import _tw, S_ST, bm_pars

def main():
    fs = 48000                  # Hz

    forward = np.loadtxt('forward.txt') * S_ST
    t = np.arange(len(forward)) / fs * 1000

    ref = np.loadtxt('uIHC_70.txt')
    sec = 99 - 70               # DSAM compatibility

    xbm = _tw.run_bm_wave(fs=48000,
                          signal=forward)

    xbm = xbm[:,sec]

    lcr4 = _tw._run_single_LCR4(fs, xbm, sec)

    uihc = _tw._run_single_ihcrp(fs, lcr4, sec)



    sec = range(400, 1600)
    t = t[sec]
    uihc = uihc[sec]
    ref = ref[sec]
    p = biggles.FramedPlot()
    p.add( biggles.Curve(t, ref, color='red', width=3) )
    p.add( biggles.Curve(t, uihc))
    p.show()




if __name__ == "__main__":
    main()
