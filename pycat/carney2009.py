# Author: Marek Rudnicki
# Time-stamp: <2009-10-09 19:27:23 marek>
#
# Description: Model of auditory periphery of: Zilany, M.S.A., Bruce,
# I.C., Nelson, P.C., and Carney, L.H. (manuscript in preparation) 2009

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import catmodel
import thorns as th

class Carney2009(object):
    def __init__(self, hsr=1, msr=1, lsr=1,
                 freq=1000, animal='cat', implnt='approx'):
        """
        hsr, msr, lsr: number of HSR/MSR/LSR fibers

        freq: CF

        animal: must be cat

        implnt: 'approx'/'acctual' implementation of the power-law
        """
        assert animal == 'cat'

        self.hsr_num = hsr
        self.msr_num = msr
        self.lsr_num = lsr

        self._implnt = implnt

        self._cohc = 1
        self._cihc = 1

        self.set_freq(freq)


    def run(self, fs, sound, times=1, output_format='spikes'):

        # TODO: implement output_format='signals'
        if output_format == 'signals':
            assert False
            assert times == 1

        if self.hsr_num > 0:
            ihc_pars = {'cohc':self._cohc, 'cihc':self._cihc}
            synapse_pars = {'nrep':self.hsr_num, 'anf_type':'hsr', 'implnt':'actual'}
            hsr_out = self._run_anf(fs, sound, times, ihc_pars, synapse_pars)
        else:
            hsr_out = None

        if self.msr_num > 0:
            ihc_pars = {'cohc':self._cohc, 'cihc':self._cihc}
            synapse_pars = {'nrep':self.msr_num, 'anf_type':'msr', 'implnt':'actual'}
            msr_out = self._run_anf(fs, sound, times, ihc_pars, synapse_pars)
        else:
            msr_out = None

        if self.lsr_num > 0:
            ihc_pars = {'cohc':self._cohc, 'cihc':self._cihc}
            synapse_pars = {'nrep':self.lsr_num, 'anf_type':'lsr', 'implnt':'actual'}
            lsr_out = self._run_anf(fs, sound, times, ihc_pars, synapse_pars)
        else:
            lsr_out = None

        return hsr_out, msr_out, lsr_out



    def _run_anf(self, fs, sound, times, ihc_pars, synapse_pars):
        anf_db = []
        for freq_idx,cf in enumerate(self._freq_map):
            vihc = catmodel.run_ihc(fs=fs, sound=sound, cf=cf, **ihc_pars)
            for run_idx in range(times):
                psth = catmodel.run_synapse(fs=fs,
                                         vihc=vihc,
                                         cf=cf,
                                         **synapse_pars);

                train = th.signal_to_spikes(fs, psth)
                train = train[0] # there is only one train per run
                anf_db.append( (freq_idx, run_idx, train) )
        anf_out = np.array(anf_db, dtype=[ ('freq', int),
                                           ('trial', int),
                                           ('spikes', np.ndarray) ])
        return anf_out

    def set_freq(self, freq):
        if isinstance(freq, int):
            freq = float(freq)
        assert (isinstance(freq, tuple) or
                isinstance(freq, float))

        if isinstance(freq, float):
            self._freq_map = [freq]
        elif isinstance(freq, tuple):
            # Based on GenerateGreenwood_CFList() from DSAM
            aA = 456.0
            k = 0.8
            a = 2.1

            freq_min, freq_max, freq_num = freq

            xmin = np.log10( freq_min / aA + k) / a
            xmax = np.log10( freq_max / aA + k) / a

            x_map = np.linspace(xmin, xmax, freq_num)
            self._freq_map = aA * ( 10**( a*x_map ) - k)
        else:
            assert False

    def get_freq_map(self):
        return self._freq_map


def main():
    ear = Carney2009(hsr=1, msr=0, lsr=0)

    fs = 100000.0
    cf = 1000
    stimdb = 50
    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = catmodel.set_dB_SPL(50, s)

    ear.set_freq( cf )

    hsr, msr, lsr = ear.run(fs, s, times=250)
    th.plot_raster(hsr['spikes'])


if __name__ == "__main__":
    main()

