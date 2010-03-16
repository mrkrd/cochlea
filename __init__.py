from __future__ import division

import numpy as np
import scipy.signal as dsp
import _tw
import bm_pars
from bm_pars import real_freq_map, S_ST, S_ED, C_eardrum


# np.fliplr are necessary for compatimility with DSAM

# Input signal should be multiplied by this factor
scaling_factor = S_ST * S_ED


def run_bm(fs, signal, mode='x', enable_LCR4=True):
    """ Simulate response of the Basilar membrane to audio siganl.

    fs: sampling frequency
    signal: audio signal sampled at 48000Hz (uPa)
    mode: 'x' return BM displacement
          'v' return BM velocity
    enable_LCR4: if False skips LCR4 compression stage

    return: BM displacement or velocity for 100 frequency sections
            (real_freq_map)

    """
    assert fs == 48000

    fs = float(fs)

    signal = np.squeeze(signal)
    assert signal.ndim == 1

    # uPa -> Pa
    signal = signal * 1e-6

    if enable_LCR4:
        # pad signal with max delay
        orig_signal_len = len(signal)
        delays = np.round(bm_pars.delay_time * fs)
        signal = np.append(signal,
                           np.zeros(delays.max()))

    _tw.bm_init(fs=fs,
                Ls=bm_pars.Ls,
                Rs=bm_pars.Rs,
                Ct=bm_pars.Ct,
                Rbm=bm_pars.Rbm,
                Cbm=bm_pars.Cbm,
                Lbm=bm_pars.Lbm,
                Rh=bm_pars.Rh, # helicotrema, end of BM
                Lh=bm_pars.Lh) # end of BM

    xBM = _tw.bm_wave(signal=signal,
                      ampl_corr=bm_pars.ampl_corr,
                      Abm=bm_pars.Abm,
                      Cbm=bm_pars.Cbm)


    if enable_LCR4:

        _tw.LCR4_init(fs=fs,
                      freq_map=bm_pars.freq_map,
                      Qmin=bm_pars.Qmin,
                      SAT1=bm_pars.SAT1,
                      SAT4=bm_pars.SAT4)


        xBM = _tw.LCR4(xBM=xBM,
                       Qmax=bm_pars.Qmax,
                       Qmin=bm_pars.Qmin)


        # Compensate for LCR4 delays
        i = np.tile(np.arange(orig_signal_len), 100).reshape( (100, orig_signal_len) ).T
        i = i + delays
        i = i.astype(int)

        j = np.arange(100)

        xBM = xBM[i,j]


    if mode == 'x':
        outBM = xBM
    elif mode == 'v':
        outBM = np.diff(xBM, axis=0) * fs

    return np.fliplr(outBM)




def run_ihcrp(fs, xBM):
    """ Run modified Shama DSAM module.

    uIHC = run_ihcrp(fs, xBM)

    fs: sampling frequency
    xBM: BM displacement

    uIHC: IHC potential
    """
    fs = float(fs)
    _tw.ihcrp_init(fs)

    xBM = np.fliplr(xBM)
    uIHC = _tw.ihcrp(xBM, bm_pars.ciliaGain)
    uIHC = np.fliplr(uIHC)

    return uIHC



def run_middle_ear_filter_orig(fs, signal):
    """ Middle ear filter designed using digital wave techinique. """
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



def run_middle_ear_filter(fs, signal):
    """ Middle ear filter model. """
    # Digital wave filter coefficients
    R2_ME = 1. / (2. * fs * C_eardrum)
    R1 = 1. / (2. * np.pi * C_eardrum * 1e3)
    g1_ME = (R2_ME - R1) / (R2_ME + R1)
    Z_ME=0

    # Standard filter coefficients
    b = [(-1-g1_ME), (1+g1_ME)]
    a = [R2_ME, R2_ME*g1_ME]

    return dsp.lfilter(b, a, signal)


def run_outer_ear_filter(fs, signal):
    assert fs == 48000

    a = bm_pars.outer_ear_a
    b = bm_pars.outer_ear_b

    return dsp.lfilter(b, a, signal)
