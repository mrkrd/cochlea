import numpy as np
import scipy.signal as dsp
import _bm
import bm_pars
from bm_pars import real_freq_map, S_ST, S_ED, C_eardrum


# np.fliplr are necessary for compatimility with DSAM

def run_bm(fs, signal, mode='v', with_LCR=True):
    """
    Simulate response of the Basilar membrane to audio siganl.

    fs: sampling frequency
    signal: audio signal sampled at 48000Hz
    mode: 'x' return BM displacement
          'v' return BM velocity
    with_LCR: if False skips LCR compression stage

    return: BM displacement or velocity for 100 frequency sections
            (real_freq_map)
    """
    assert fs == 48000

    fs = float(fs)

    signal = np.squeeze(signal)
    assert signal.ndim == 1

    # uPa -> Pa
    signal = signal * 1e-6

    _bm.bm_init(fs,
                bm_pars.Ls,
                bm_pars.Rs,
                bm_pars.Ct,
                bm_pars.Rbm,
                bm_pars.Cbm,
                bm_pars.Lbm,
                bm_pars.Rh, # helicotrema, end of BM
                bm_pars.Lh) # end of BM

    xBM = _bm.bm_wave(signal,
                      bm_pars.ampl_corr,
                      bm_pars.Abm,
                      bm_pars.Cbm)


    if with_LCR:

        _bm.LCR4_init(fs,
                      bm_pars.freq_map,
                      bm_pars.Qmin,
                      bm_pars.SAT1,
                      bm_pars.SAT4)


        xBM = _bm.LCR4(xBM,
                       bm_pars.Qmax,
                       bm_pars.Qmin);


    if mode == 'x':
        outBM = xBM
    elif mode == 'v':
        outBM = np.diff(xBM, axis=0) * fs

    return np.fliplr(outBM)




def run_ihcrp(fs, xBM):
    """
    Run modified Shama DSAM module.

    uIHC = run_ihcrp(fs, xBM)

    fs: sampling frequency
    xBM: BM displacement

    uIHC: IHC potential
    """
    fs = float(fs)
    _bm.ihcrp_init(fs)

    uIHC = _bm.ihcrp(np.fliplr(xBM), bm_pars.ciliaGain)

    return np.fliplr(uIHC)



def run_stapes(s):
    """
    TODO: docs
    """

    # 5.61382474984 is the abs value of the orignal filter's
    # transferfunction at 1kHz
    return s * S_ST * S_ED * 5.61382474984


def run_mid_ear_filter_orig(fs, signal):
    """
    Middle ear filter designed using digital wave techinique.
    """
    R2_ME = 1. / (2. * fs * C_eardrum)
    R1 = 1. / (2. * np.pi * C_eardrum * 1e3)

    g1_ME = (R2_ME - R1) / (R2_ME + R1)
    Z_ME=0

    out = np.zeros(len(signal))

    for i,samp in enumerate(signal):
        b2 = samp + g1_ME * (samp - Z_ME)
        Z_ME_b = Z_ME
        Z_ME = b2
        out[i] = ((Z_ME_b - b2) / R2_ME)

    return out


def run_mid_ear_filter(fs, signal):
    """
    Middle ear filter from Werner's model.
    """
    # Digital wave filter coefficients
    R2_ME = 1. / (2. * fs * C_eardrum)
    R1 = 1. / (2. * np.pi * C_eardrum * 1e3)
    g1_ME = (R2_ME - R1) / (R2_ME + R1)
    Z_ME=0

    # Standard filter coefficients
    b = [(-1-g1_ME), (1+g1_ME)]
    a = [R2_ME, R2_ME*g1_ME]

    return dsp.lfilter(b, a, signal)
