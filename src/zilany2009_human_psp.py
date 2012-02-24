# Description: Model of auditory periphery of: Zilany, M.S.A., Bruce,
# I.C., Nelson, P.C., and Carney, L.H. (manuscript in preparation) 2009

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import scipy.signal as dsp
import os

import _pycat

import thorns as th


class Zilany2009_Human_PSP(object):
    name = 'Zilany2009_Human_PSP'

    def __init__(self, anf_num=(1,1,1), cf=1000,
                 powerlaw_implnt='approx', with_ffGn=False,
                 with_me=True):
        """Auditory periphery model from Zilany et al. (2009) with
        human middle ear filter

        anf_num: (hsr_num, msr_num, lsr_num)
        cf: CF
        powerlaw_implnt: 'approx' or 'actual' implementation of the power-law
        with_ffGn: enable/disable Gausian noise
        with_me: enable/disable middle ear filter (used for fitting a new ME filter)

        """
        assert not set(anf_num) - set((1,0)), "anf_num can have only 0 and 1 values"


        self._hsr_num = anf_num[0]
        self._msr_num = anf_num[1]
        self._lsr_num = anf_num[2]


        self._powerlaw_implnt = powerlaw_implnt
        self._with_ffGn = with_ffGn

        self._cohc = 1
        self._cihc = 1

        self._with_me = with_me

        self.set_freq(cf)


    def run(self, sound, fs, seed):
        """ Run the model.

        fs: sampling frequency of the signal; model is run at the same frequency
        sound: input signal

        """
        np.random.seed(seed)

        assert np.max(sound) < 1000, "Signal should be given in Pa"


        # Run Outer/Middle Ear filter
        if self._with_me:
            sound = run_me_filter_for_zilany2009(sound, fs)


        psp = {'hsr':[], 'msr':[], 'lsr':[]}
        for cf in self._freq_map:

            # Run IHC model
            vihc = _pycat.run_ihc(signal=sound, cf=cf, fs=fs,
                                  cohc=self._cohc, cihc=self._cihc)


            # Run HSR synapse
            if self._hsr_num:
                p = self._run_anf(fs, cf, vihc,
                                  anf_type='hsr',
                                  anf_num=self._hsr_num)
                psp['hsr'].append( p )

            # Run MSR synapse
            if self._msr_num:
                p = self._run_anf(fs, cf, vihc,
                                  anf_type='msr',
                                  anf_num=self._msr_num)
                psp['msr'].append( p )

            # Run LSR synapse
            if self._lsr_num:
                p = self._run_anf(fs, cf, vihc,
                                  anf_type='lsr',
                                  anf_num=self._lsr_num)
                psp['lsr'].append( p )


        for typ in psp:
            psp[typ] = np.array(psp[typ]).T

        return psp


    def _run_anf(self, fs, cf, vihc, anf_type, anf_num):

        duration = len(vihc) / fs # [s]

        synout = _pycat.run_synapse(fs=fs, vihc=vihc, cf=cf,
                                    anf_type=anf_type,
                                    powerlaw_implnt=self._powerlaw_implnt,
                                    with_ffGn=self._with_ffGn)

        return synout


    def set_freq(self, cf):
        """ Set signle or range of CF for the model."""
        if isinstance(cf, int):
            cf = float(cf)
        assert (isinstance(cf, tuple) or
                isinstance(cf, float))

        if isinstance(cf, float):
            self._freq_map = [cf]
        elif isinstance(cf, tuple):
            # Based on GenerateGreenwood_CFList() from DSAM
            # Liberman (1982)
            # Pars from GetGreenwoodPars_CFList()
            aA = 165.4
            k = 0.88
            a = 2.1

            freq_min, freq_max, freq_num = cf

            xmin = np.log10( freq_min / aA + k) / a
            xmax = np.log10( freq_max / aA + k) / a

            x_map = np.linspace(xmin, xmax, freq_num)
            self._freq_map = aA * ( 10**( a*x_map ) - k)

    def get_freq_map(self):
        return self._freq_map




def run_me_filter_for_zilany2009(signal, fs):
    assert fs > 40e3

    signal_fft = np.fft.fft(signal)
    freqs = np.fft.fftfreq(len(signal), d=1/fs)


    me_filter_response = np.loadtxt(data_dir("me_filter_response.txt"))
    me_filter_freqs = np.loadtxt(data_dir("me_filter_freqs.txt"))
    fmin = me_filter_freqs.min()
    fmax = me_filter_freqs.max()


    # Convert dB to amplitudae ratio
    response_ratio = 10 ** (me_filter_response / 20)

    # Interpolate the filter to fit signal's FFT
    band_len = len(signal_fft[(freqs>=fmin) & (freqs<=fmax)])
    interp_ratio = dsp.resample(response_ratio, band_len)
    interp_ratio = np.concatenate( (interp_ratio, interp_ratio[::-1]) )


    # Apply the filter
    band = (np.abs(freqs) >= fmin) & (np.abs(freqs) <= fmax)
    signal_fft[band] *= interp_ratio

    signal_fft[ np.logical_not(band) ] = 0



    filtered = np.fft.ifft( signal_fft )

    return np.real(filtered)


def data_dir(par_file):
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'data', par_file)



def main():
    import thorns as th
    import thorns.waves as wv

    import matplotlib.pyplot as plt

    fs = 100e3
    cf = (80, 8000, 100)
    stimdb = 10

    ear = Zilany2009_Human_PSP((1,1,1), cf=cf,
                               powerlaw_implnt='approx',
                               with_ffGn=False
                           )

    s = wv.generate_ramped_tone(fs,
                                freq=1000,
                                tone_duration=50e-3,
                                ramp_duration=2.5e-3,
                                pad_duration=20e-3,
                                dbspl=stimdb)

    psp = ear.run(s, fs, None)

    plt.imshow(psp['hsr'], aspect='auto')
    plt.show()



if __name__ == "__main__":
    main()

