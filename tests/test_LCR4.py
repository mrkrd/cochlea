from __future__ import division

import numpy as np
import biggles

from traveling_waves import _tw, S_ST, bm_pars

def main():
    fs = 48000                  # Hz

    forward = np.loadtxt('forward.txt') * S_ST
    t = np.arange(len(forward)) / fs * 1000

    ref = np.loadtxt('LCR4_70.txt')


    _tw.bm_init(fs,
                bm_pars.Ls,
                bm_pars.Rs,
                bm_pars.Ct,
                bm_pars.Rbm,
                bm_pars.Cbm,
                bm_pars.Lbm,
                bm_pars.Rh, # helicotrema, end of BM
                bm_pars.Lh) # end of BM
    xBM = _tw.bm_wave(forward,
                      bm_pars.ampl_corr,
                      bm_pars.Abm,
                      bm_pars.Cbm)
    _tw.LCR4_init(fs,
                  bm_pars.freq_map,
                  bm_pars.Qmin,
                  bm_pars.SAT1,
                  bm_pars.SAT4)
    LCR4 = _tw.LCR4(xBM,
                    bm_pars.Qmax,
                    bm_pars.Qmin)


    LCR4 = LCR4[:,70]


    sec = range(400, 2000)
    t = t[sec]
    LCR4 = LCR4[sec]
    ref = ref[sec]
    p = biggles.FramedPlot()
    p.add( biggles.Curve(t, ref, color='red', width=3) )
    p.add( biggles.Curve(t, LCR4))
    p.show()




if __name__ == "__main__":
    main()
