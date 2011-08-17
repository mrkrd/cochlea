#!/usr/bin/env python

"""Holmberg, M. (2007). Speech Encoding in the Human Auditory
Periphery: Modeling and Quantitative Assessment by Means of Automatic
Speech Recognition. PhD thesis, Technical University Darmstadt.

"""


from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

from auditory_periphery import par_dir, AuditoryPeriphery
import traveling_waves as tw
import dsam


class Holmberg2007(AuditoryPeriphery):
    name = "Holmberg2007"

    def __init__(self, anf_num=(1,1,1),
                 cf=None,
                 accumulate=False,
                 approx_freq=False):
        """ Auditory periphery model from Marcus Holmberg (2007)

        anf_num: (hsr_num, msr_num, lsr_num)
        cf: CF (if None all 100 channels are computed)
        accumulate: if True, spikes for all fibers are calculated at once.
        approx_freq: if True, calculate channel that is the closest to the given `cf'
                     No assertion error will be generated.

        """

        self.approx_freq = approx_freq
        self.set_freq(cf)

        self._hsr_num = anf_num[0]
        self._msr_num = anf_num[1]
        self._lsr_num = anf_num[2]

        self._accumulate = accumulate

        # Outer/middle ear filter
        # Stapes velocity [Pa -> m/s]
        # BM module is in run()
        # IHC receptor potential


        if self._hsr_num > 0:
            self.ihc_hsr_module = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr_module.read_pars(par_dir("ihc_hsr_Sumner2002.par"))

            self.sg_hsr_module = dsam.EarModule("An_SG_Carney")
            self.sg_hsr_module.read_pars(par_dir("anf_carney.par"))
            dsam.connect(self.ihc_hsr_module, self.sg_hsr_module)

        if self._msr_num > 0:
            self.ihc_msr_module = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr_module.read_pars(par_dir("ihc_msr_Sumner2002.par"))

            self.sg_msr_module = dsam.EarModule("An_SG_Carney")
            self.sg_msr_module.read_pars(par_dir("anf_carney.par"))
            dsam.connect(self.ihc_msr_module, self.sg_msr_module)


        if self._lsr_num > 0:
            self.ihc_lsr_module = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr_module.read_pars(par_dir("ihc_lsr_Sumner2002.par"))

            self.sg_lsr_module = dsam.EarModule("An_SG_Carney")
            self.sg_lsr_module.read_pars(par_dir("anf_carney.par"))
            dsam.connect(self.ihc_lsr_module, self.sg_lsr_module)


    def set_freq(self, cf):

        if isinstance(cf, int):
            cf = float(cf)

        if isinstance(cf, float):
            if self.approx_freq:
                self._freq_idx = tw.find_closest_freq_idx_in_map(cf)
            else:
                real_freq_map = tw.bm_pars.real_freq_map
                assert cf in real_freq_map
                self._freq_idx = int(np.where(real_freq_map == cf)[0])
        elif cf is None:
            self._freq_idx = None
        else:
            assert False, "CF must be a real number or None"



    def get_freq_map(self):
        """ Returns frequency map of the model. """
        if isinstance(self._freq_idx, int):
            freq_map = [ tw.bm_pars.real_freq_map[self._freq_idx] ]
        else:
            freq_map = tw.bm_pars.real_freq_map

        return freq_map


    def run(self, fs, sound):
        """ Run auditory periphery model.

        sound: audio signal
        fs: sampling frequency (48000 Hz)

        """
        assert fs == 48000

        ### Outer ear filter
        sound = tw.run_outer_ear_filter(fs, sound)

        ### Middle ear filter
        sound = tw.run_middle_ear_filter(fs, sound)

        ### Scaling
        sound = sound * tw.scaling_factor

        ### Basilar membrane
        xbm = tw.run_bm_wave(fs, sound)

        if self._freq_idx is not None:
            xbm = xbm[:,self._freq_idx]

        ### Amplification
        LCR4 = tw.run_LCR4(fs, xbm, self._freq_idx)

        ### IHCRP
        ihcrp = tw.run_ihcrp(fs, LCR4, self._freq_idx)


        trains = []
        if self._hsr_num > 0:
            self.ihc_hsr_module.run(fs, ihcrp)
            tr = self._run_anf('hsr', self.sg_hsr_module,
                               fs, self._hsr_num, self._accumulate)
            trains.extend(tr)

        if self._msr_num > 0:
            self.ihc_msr_module.run(fs, ihcrp)
            tr = self._run_anf('msr', self.sg_msr_module,
                               fs, self._msr_num, self._accumulate)
            trains.extend(tr)

        if self._lsr_num > 0:
            self.ihc_lsr_module.run(fs, ihcrp)
            tr = self._run_anf('lsr', self.sg_lsr_module,
                               fs, self._lsr_num, self._accumulate)
            trains.extend(tr)


        spike_trains = np.rec.array(trains, dtype=self._anf_dtype)
        return spike_trains



    @classmethod
    def plot_rate_intensity(cls):
        ear = Holmberg2007()
        print ear


def main():
    import thorns as th
    import thorns.waves as wv

    fs = 48000
    cf = tw.find_closest_freq_in_map(1000)
    print "CF:", cf
    stimdb = 60

    ear = Holmberg2007((250,0,0), cf=cf)


    s = wv.generate_ramped_tone(fs,
                                freq=cf,
                                tone_duration=50,
                                ramp_duration=2.5,
                                pad_duration=20,
                                dbspl=stimdb)

    anf = ear.run(fs, s)

    th.plot.raster(anf).show()
    th.plot.psth(anf, bin_size=1).show()






if __name__ == "__main__":
    main()
