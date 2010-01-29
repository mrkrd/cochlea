from __future__ import division

import numpy as np

from auditory_periphery import par_dir, AuditoryPeriphery
import traveling_waves as tw
import dsam


class Holmberg2008(AuditoryPeriphery):
    def __init__(self, anf_sum=(100, 100, 100),
                 freq=None, animal='human', sg_type='carney'):
        """ Auditory periphery model from Marcus Holmberg

        anf_sum: (hsr_sum, msr_sum, lsr_sum)
        freq: CF
        animal: only 'human' implemented

        """

        assert animal == 'human'

        self.set_freq(freq)

        self._hsr_sum = anf_sum[0]
        self._msr_sum = anf_sum[1]
        self._lsr_sum = anf_sum[2]
        self._animal = animal


        # Outer/middle ear filter
        self.outer_middle_ear = dsam.EarModule("Filt_MultiBPass")
        self.outer_middle_ear.read_pars(par_dir("filt_Human.par"))


        # Stapes velocity [Pa -> m/s]
        # BM module is in run()
        # IHC receptor potential


        if self._hsr_sum != 0:
            self.ihc_hsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr.read_pars(par_dir("ihc_hsr_Meddis2002.par"))

            self.anf_hsr = self._generate_anf(sg_type, self._hsr_sum)
            dsam.connect(self.ihc_hsr, self.anf_hsr)

        if self._msr_sum != 0:
            self.ihc_msr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr.read_pars(par_dir("ihc_msr_Meddis2002.par"))

            self.anf_msr = self._generate_anf(sg_type, self._msr_sum)
            dsam.connect(self.ihc_msr, self.anf_msr)


        if self._lsr_sum != 0:
            self.ihc_lsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr.read_pars(par_dir("ihc_lsr_Meddis2002.par"))

            self.anf_lsr = self._generate_anf(sg_type, self._lsr_sum)
            dsam.connect(self.ihc_hsr, self.anf_lsr)


    def set_freq(self, freq):

        if isinstance(freq, int):
            freq = float(freq)

        # Only real numbers please.
        assert (isinstance(freq, float) or
                freq == None)

        if isinstance(freq, float):
            real_freq_map = tw.bm_pars.real_freq_map
            assert freq in real_freq_map
            self._freq_idx = int(np.where(real_freq_map == freq)[0])
        elif freq == None:
            self._freq_idx = None


    def get_freq_map(self):
        """ Returns frequency map of the model. """
        if isinstance(self._freq_idx, int):
            freq_map = [ tw.bm_pars.real_freq_map[self._freq_idx] ]
        else:
            freq_map = tw.bm_pars.real_freq_map

        return freq_map


    def run(self, fs, sound, times=1):
        """ Run auditory periphery model.

        sound: audio signal
        fs: sampling frequency
        times: how many many trials

        """
        # TODO: filter with original filters from the model
        self.outer_middle_ear.run(fs, sound)

        s = self.outer_middle_ear.get_signal()

        ### Stapes
        stapes_velocity = tw.run_stapes(s)

        ### Basilar membrane
        xBM = tw.run_bm(fs, stapes_velocity, mode='x')

        ### IHCRP
        ihcrp = tw.run_ihcrp(fs, xBM)
        if self._freq_idx != None:
            ihcrp = ihcrp[:,self._freq_idx]


        if self._hsr_sum > 0:
            self.ihc_hsr.run(fs, ihcrp)
            hsr_trains = self._run_anf(self.anf_hsr, fs, times)
        else:
            hsr_trains = None

        if self._msr_sum > 0:
            self.ihc_msr.run(fs, ihcrp)
            msr_trains = self._run_anf(self.anf_msr, fs, times)
        else:
            msr_trains = None

        if self._lsr_sum > 0:
            self.ihc_lsr.run(fs, ihcrp)
            lsr_trains = self._run_anf(self.anf_lsr, fs, times)
        else:
            lsr_trains = None

        return hsr_trains, msr_trains, lsr_trains



def main():
    import thorns as th

    fs = 48000
    cf = tw.real_freq_map[38]

    ear = Holmberg2008((1,0,0), freq=cf)

    stimdb = 50

    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = dsam.set_dbspl(stimdb, s)
    z = np.zeros( np.ceil(len(t)/2) )
    s = np.concatenate( (z, s, z) )

    ear.set_freq( cf )

    hsr, msr, lsr = ear.run(fs, s, times=250)
    th.plot_raster(hsr['spikes'])
    th.plot_psth(hsr['spikes'])


if __name__ == "__main__":
    main()
