# Author: Marek Rudnicki
# Time-stamp: <2009-12-18 00:32:03 marek>
#
# Description: Model of auditory periphery of: Zilany, M.S.A., Bruce,
# I.C., Nelson, P.C., and Carney, L.H. (manuscript in preparation) 2009

from __future__ import division

import numpy as np

import catmodel
import thorns as th


class Zilany2009(object):
    def __init__(self, anf_sum=(1,1,1), freq=1000,
                 animal='cat', powerlaw_implnt='actual'):
        """Auditory periphery model of a cat (Zilany et al. 2009)

        anf_sum: (hsr_sum, msr_sum, lsr_sum)
        freq: CF
        animal: must be 'cat'
        powerlaw_implnt: 'approx' or 'acctual' implementation of the power-law

        """
        assert animal == 'cat'

        self._hsr_sum = anf_sum[0]
        self._msr_sum = anf_sum[1]
        self._lsr_sum = anf_sum[2]

        self._powerlaw_implnt = powerlaw_implnt

        self._cohc = 1
        self._cihc = 1

        self.set_freq(freq)


    def run(self, fs, sound, times=1):
        """ Run the model.

        fs: sampling frequency of the signal; model is run at the same frequency
        sound: input signal
        times: number of tirals / number of output trains

        """
        # TODO: implement storing the spikes in a file and reloading them automaticly

        hsr_trains = []
        msr_trains = []
        lsr_trains = []
        for cf in self._freq_map:
            # Run IHC model
            vihc = catmodel.run_ihc(fs=fs, sound=sound, cf=cf, cohc=self._cohc, cihc=self._cihc)

            # Run HSR synapse
            if self._hsr_sum > 0:
                synapse_pars = {'anf_type':'hsr',
                                'anf_sum':self._hsr_sum,
                                'powerlaw_implnt':self._powerlaw_implnt}
                hsr_trains.extend( self._run_anf(fs, cf, vihc, times, synapse_pars) )
            else:
                hsr_trains = None

            # Run MSR synapse
            if self._msr_sum > 0:
                synapse_pars = {'anf_type':'msr',
                                'anf_sum':self._msr_sum,
                                'powerlaw_implnt':self._powerlaw_implnt}
                msr_trains.extend( self._run_anf(fs, cf, vihc, times, synapse_pars) )
            else:
                msr_trains = None

            # Run LSR synapse
            if self._lsr_sum > 0:
                synapse_pars = {'anf_type':'lsr',
                                'anf_sum':self._lsr_sum,
                                'powerlaw_implnt':self._powerlaw_implnt}
                lsr_trains.extend( self._run_anf(fs, cf, vihc, times, synapse_pars) )
            else:
                lsr_trains = None


        train_type = [ ('freq', float), ('trial', int), ('spikes', np.ndarray) ]
        hsr_trains = np.array(hsr_trains, dtype=train_type)
        msr_trains = np.array(msr_trains, dtype=train_type)
        lsr_trains = np.array(lsr_trains, dtype=train_type)

        return hsr_trains, msr_trains, lsr_trains



    def _run_anf(self, fs, cf, vihc, times, synapse_pars):

        anf_trains = []
        for run_idx in range(times):
            psth = catmodel.run_synapse(fs=fs,
                                        vihc=vihc,
                                        cf=cf,
                                        **synapse_pars);
            train = th.signal_to_spikes(fs, psth)
            train = train[0] # there is only one train per run
            anf_trains.append( (cf, run_idx, train) )

        return anf_trains


    def set_freq(self, freq):
        """ Set signle or range of CF for the model."""
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

    def get_freq_map(self):
        return self._freq_map



def main():
    ear = Zilany2009((1,0,0), freq=5000, powerlaw_implnt='approx')

    fs = 100000.0
    cf = 1000
    stimdb = 50

    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = catmodel.set_dbspl(stimdb, s)
    z = np.zeros( np.ceil(len(t)/2) )
    s = np.concatenate( (z, s, z) )

    ear.set_freq( cf )

    hsr, msr, lsr = ear.run(fs, s, times=250)
    th.plot_raster(hsr['spikes'])
    th.plot_psth(hsr['spikes'])


if __name__ == "__main__":
    main()

