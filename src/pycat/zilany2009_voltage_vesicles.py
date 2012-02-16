# Description: Model of auditory periphery of: Zilany, M.S.A., Bruce,
# I.C., Nelson, P.C., and Carney, L.H. (manuscript in preparation) 2009

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import thorns as th
import _pycat


class Zilany2009_Voltage_Vesicles(object):
    def __init__(self, anf_num=(1,1,1),
                 powerlaw_implnt='approx', with_ffGn=False):
        """ Auditory periphery model of a cat (Zilany et al. 2009)

        anf_num: (hsr_num, msr_num, lsr_num)
        powerlaw_implnt: 'approx' or 'actual' implementation of the power-law

        """
        self.name = 'Zilany2009'

        self._hsr_num = anf_num[0]
        self._msr_num = anf_num[1]
        self._lsr_num = anf_num[2]

        self._powerlaw_implnt = powerlaw_implnt
        self._with_ffGn = with_ffGn

        self._cohc = 1
        self._cihc = 1

        self._train_type = [('typ', 'S3'),
                            ('cf', float),
                            ('spikes', np.ndarray)]


    def run(self, fs, voltage):
        """ Run the model.

        fs: sampling frequency of the signal; model is run at the same frequency
        voltage: IHC receptor potential (V)

        """
        voltage = voltage + 55e-3
        trains = []

        # Run HSR synapse
        if self._hsr_num > 0:
            tr = self._run_anf(fs, voltage,
                               anf_type='hsr',
                               anf_num=self._hsr_num)
            trains.extend(tr)

        # Run MSR synapse
        if self._msr_num > 0:
            tr = self._run_anf(fs, voltage,
                               anf_type='msr',
                               anf_num=self._msr_num)
            trains.extend(tr)

        # Run LSR synapse
        if self._lsr_num > 0:
            tr = self._run_anf(fs, voltage,
                               anf_type='lsr',
                               anf_num=self._lsr_num)
            trains.extend(tr)

        trains = np.array(trains, dtype=self._train_type)

        return trains



    def _run_anf(self, fs, vihc, anf_type, anf_num):

        synout = None
        anf_trains = []
        for anf_id in range(anf_num):
            if (synout is None) or (self._with_ffGn):
                synout = _pycat.run_synapse(fs=fs, vihc=vihc, cf=1000,
                                            anf_type=anf_type,
                                            powerlaw_implnt=self._powerlaw_implnt,
                                            with_ffGn=self._with_ffGn)

            vesicles = _pycat.sample_rate(fs, synout)
            vesicles = th.signal_to_spikes(fs, vesicles)[0]

            anf_trains.append( (anf_type, -1, vesicles) )

        return anf_trains



def main():
    fs = 100000
    cf = 100

    ear = Zilany2009_Voltage_Vesicles((1000,0,0),
                                      powerlaw_implnt='approx',
                                      with_ffGn=False)

    n = np.ceil(0.1*fs)
    v = np.ones(n) * -89e-3                 # V
    v[np.ceil(n/3):np.ceil(n*2/3)] = -29e-3 # V

    anf = ear.run(fs, v)

    th.plot_raster(anf['spikes']).show()

    p = th.plot_psth(anf[anf['typ']=='hsr']['spikes'], color='black')
    p.show()


if __name__ == "__main__":
    main()

