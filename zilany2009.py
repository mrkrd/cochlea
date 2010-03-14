# Author: Marek Rudnicki
# Time-stamp: <2010-03-14 15:11:15 marek>
#
# Description: Model of auditory periphery of: Zilany, M.S.A., Bruce,
# I.C., Nelson, P.C., and Carney, L.H. (manuscript in preparation) 2009

from __future__ import division

import numpy as np

import catmodel
import thorns as th


class Zilany2009(object):
    def __init__(self, anf_num=(1,1,1), cf=1000,
                 powerlaw_implnt='actual', accumulate=False, verbose=False):
        """ Auditory periphery model of a cat (Zilany et al. 2009)

        anf_num: (hsr_num, msr_num, lsr_num)
        cf: CF
        powerlaw_implnt: 'approx' or 'acctual' implementation of the power-law
        accumulate: if True, then spike trains of each type are concatenated

        """
        self.name = 'Zilany2009'

        assert accumulate == False

        self._hsr_num = anf_num[0]
        self._msr_num = anf_num[1]
        self._lsr_num = anf_num[2]

        self._powerlaw_implnt = powerlaw_implnt

        self._cohc = 1
        self._cihc = 1

        self._train_type = [('typ', 'S3'),
                            ('cf', float),
                            ('id', int),
                            ('spikes', np.ndarray)]

        self.set_freq(cf)

        self._verbose = verbose

    def run(self, fs, sound):
        """ Run the model.

        fs: sampling frequency of the signal; model is run at the same frequency
        sound: input signal

        """
        # TODO: implement storing of spikes in a file/db and reloading them as needed

        trains = []
        for cf in self._freq_map:
            # Run IHC model
            vihc = catmodel.run_ihc(fs=fs, sound=sound, cf=cf,
                                    cohc=self._cohc, cihc=self._cihc,
                                    verbose=self._verbose)

            # Run HSR synapse
            if self._hsr_num > 0:
                tr = self._run_anf(fs, cf, vihc,
                                   anf_type='hsr',
                                   anf_num=self._hsr_num)
                trains.extend(tr)

            # Run MSR synapse
            if self._msr_num > 0:
                tr = self._run_anf(fs, cf, vihc,
                                   anf_type='msr',
                                   anf_num=self._msr_num)
                trains.extend(tr)

            # Run LSR synapse
            if self._lsr_num > 0:
                tr = self._run_anf(fs, cf, vihc,
                                   anf_type='lsr',
                                   anf_num=self._lsr_num)
                trains.extend(tr)

        trains = np.rec.array(trains, dtype=self._train_type)

        return trains



    def _run_anf(self, fs, cf, vihc, anf_type, anf_num):

        anf_trains = []
        for anf_id in range(anf_num):
            psth = catmodel.run_synapse(fs=fs, vihc=vihc, cf=cf,
                                        anf_type=anf_type,
                                        powerlaw_implnt=self._powerlaw_implnt,
                                        verbose=self._verbose)
            train = th.signal_to_spikes(fs, psth)
            train = train[0] # there is only one train per run
            anf_trains.append( (anf_type, cf, anf_id, train) )

        return anf_trains


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
            aA = 456.0
            k = 0.8
            a = 2.1

            freq_min, freq_max, freq_num = cf

            xmin = np.log10( freq_min / aA + k) / a
            xmax = np.log10( freq_max / aA + k) / a

            x_map = np.linspace(xmin, xmax, freq_num)
            self._freq_map = aA * ( 10**( a*x_map ) - k)

    def get_freq_map(self):
        return self._freq_map



def main():
    fs = 100000.0
    cf = 1000
    stimdb = 80

    ear = Zilany2009((100,100,100), cf=cf, powerlaw_implnt='approx',
                     verbose=False)

    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = catmodel.set_dbspl(stimdb, s)
    z = np.zeros( np.ceil(len(t)/2) )
    s = np.concatenate( (z, s, z) )

    anf = ear.run(fs, s)

    th.plot_raster(anf.spikes).show()

    p = th.plot_psth(anf[anf.typ=='hsr']['spikes'], color='black')
    th.plot_psth(anf[anf.typ=='msr']['spikes'], color='red', plot=p)
    th.plot_psth(anf[anf.typ=='lsr']['spikes'], color='blue', plot=p)
    p.show()


if __name__ == "__main__":
    main()

