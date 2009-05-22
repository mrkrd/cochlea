import numpy as np
import scipy.signal as dsp
import _bm

bm_pars = np.load('bm_pars.npz')
freq_map = bm_pars['freq_map']

S_ST = 3.1e-6                          # Area of stapes [m^2]
S_ED = 55.e-6;                         # Ear drum area in m^2
C_eardrum = (0.7e-9/20.e-3/S_ED);      # Ear drum compliance nm/dbspl/m^2


def run_bm(fs, signal, mode='v', with_LCR=True):
    """
    Simulate response of the Basilar membrane to audio siganl.

    fs: sampling frequency
    signal: audio signal sampled at 48000Hz
    mode: 'x' return BM displacement
          'v' return BM velocity
    with_LCR: if False skips LCR compression stage

    return: BM displacement or velocity for 100 frequency sections
            (freq_map)
    """
    assert fs == 48000

    # Pa -> uPa
    signal = signal * 1e-6

    _bm.bm_init(48000,
                bm_pars['Ls'],
                bm_pars['Rs'],
                bm_pars['Ct'],
                bm_pars['Rbm'],
                bm_pars['Cbm'],
                bm_pars['Lbm'],
                float(bm_pars['Rh']), # helicotrema, end of BM
                float(bm_pars['Lh'])) # end of BM

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
