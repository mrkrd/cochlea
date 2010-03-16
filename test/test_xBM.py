from __future__ import division

import numpy as np
import biggles

from traveling_waves import run_middle_ear_filter, _tw, S_ST, bm_pars

def main():
    fs = 48000                  # Hz

    forward = np.loadtxt('forward.txt') * S_ST
    t = np.arange(len(forward)) / fs * 1000

    ref = np.loadtxt('xBM_70.txt')


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
    xBM = xBM[:,70]


    sec = range(400, 1000)
    t = t[sec]
    xBM = xBM[sec]
    ref = ref[sec]
    p = biggles.FramedPlot()
    p.add( biggles.Curve(t, ref, color='red', width=3))
    p.add( biggles.Curve(t, xBM))
    p.show()




if __name__ == "__main__":
    main()
