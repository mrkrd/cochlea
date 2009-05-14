import numpy as np
import _bm

bm_pars = np.load('bai_bm/bm_pars.npz')

fs = 48000.0

freq_map = bm_pars['freq_map']


def run_bm(signal, mode='x', with_LCR=True):

    _bm.bm_init(48000,
                bm_pars['Ls'],
                bm_pars['Rs'],
                bm_pars['Ct'],
                bm_pars['Rbm'],
                bm_pars['Cbm'],
                bm_pars['Lbm'],
                float(bm_pars['Rh']),
                float(bm_pars['Lh']))

    xBM = _bm.bm_wave(signal,
                      bm_pars['ampl_corr'],
                      bm_pars['Abm'],
                      bm_pars['Cbm'])


    if with_LCR:

        _bm.LCR4_init(fs,
                      bm_pars['freq_map'],
                      bm_pars['Qmin'],
                      bm_pars['SAT1'],
                      bm_pars['SAT4'])


        xBM = _bm.LCR4(xBM,
                       bm_pars['Qmax'],
                       bm_pars['Qmin']);


    if mode == 'x':
        outBM = xBM
    elif mode == 'v':
        outBM = np.diff(xBM, axis=0) * fs


    return np.fliplr(outBM)
