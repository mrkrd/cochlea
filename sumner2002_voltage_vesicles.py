#!/usr/bin/env python

"""Sumner, C. J., Lopez-Poveda, E. A., O'Mard, L. P., and Meddis,
R. (2002). A revised model of the inner-hair cell and auditory-nerve
complex. The Journal of the Acoustical Society of America,
111(5):2178-2188.

Input is the IHC depolarization (V) and output are vesicle timings.

"""

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import thorns as th
import dsam

from auditory_periphery import AuditoryPeriphery, par_dir

class Sumner2002_Voltage_Vesicles(AuditoryPeriphery):
    def __init__(self, anf_num=(1,1,1)):
        """ Auditory periphery model from Sumner et al. (2002) that
        outputs vesicle timings instead of spikes.  Input is the IHC
        depolarization in volts.

        anf_num: (hsr_num, msr_num, lsr_num)

        """
        self.name = 'Sumner2002_Voltage_Vesicles'

        self._hsr_num = anf_num[0]
        self._msr_num = anf_num[1]
        self._lsr_num = anf_num[2]



        if self._hsr_num > 0:
            self.ihc_hsr_module = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr_module.read_pars(par_dir("ihc_hsr_Sumner2002.par"))
            self.ihc_hsr_module.set_par('OP_MODE', 'spike')

        if self._msr_num > 0:
            self.ihc_msr_module = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr_module.read_pars(par_dir("ihc_msr_Sumner2002.par"))
            self.ihc_msr_module.set_par('OP_MODE', 'spike')

        if self._lsr_num > 0:
            self.ihc_lsr_module = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr_module.read_pars(par_dir("ihc_lsr_Sumner2002.par"))
            self.ihc_lsr_module.set_par('OP_MODE', 'spike')



    def run(self, voltage, fs):
        """ Run auditory periphery model.

        fs: sampling frequency
        voltage: IHC receptor potential

        """
        trains = []
        if self._hsr_num > 0:
            tr = self._run_ihc('hsr', self.ihc_hsr_module,
                               fs, voltage, self._hsr_num)
            trains.extend(tr)

        if self._msr_num > 0:
            tr = self._run_ihc('msr', self.ihc_msr_module,
                               fs, voltage, self._msr_num)
            trains.extend(tr)

        if self._lsr_num > 0:
            tr = self._run_ihc('lsr', self.ihc_lsr_module,
                               fs, voltage, self._lsr_num)
            trains.extend(tr)


        spike_trains = np.rec.array(trains, dtype=self._anf_dtype)
        return spike_trains




    def _run_ihc(self, anf_type, ihc_module, fs, voltage, anf_num):
        """ Skip spike generator several times and format the output. """

        anf_trains = []
        for anf_id in range(anf_num):
            ihc_module.run(fs, voltage)
            vesicle_signal = ihc_module.get_signal()
            vesicles = th.signal_to_spikes(fs, vesicle_signal)

            for train in vesicles:
                anf_trains.append( (anf_type, -1, train) )

        return anf_trains



def main():
    fs = 100000                 # Hz
    tmax = 1000                 # ms

    ear = Sumner2002_Voltage_Vesicles((250,0,0))

    v = -0.089 * np.ones( np.ceil(tmax*fs/1000) )
    v[np.ceil(10*fs/1000):np.ceil(990*fs/1000)] = -0.029

    anf = ear.run(v, fs)

    p = th.plot_raster(anf['spikes'])
    p.show()
    p = th.plot_psth(anf['spikes'])
    p.show()




if __name__ == "__main__":
    main()

