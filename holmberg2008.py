from __future__ import division

import numpy as np

from auditory_periphery import par_dir, AuditoryPeriphery
import traveling_waves as tw
import dsam


class Holmberg2008(AuditoryPeriphery):
    def __init__(self, anf_num=(100, 100, 100),
                 cf=None, sg_type='carney'):
        """ Auditory periphery model from Marcus Holmberg

        anf_num: (hsr_num, msr_num, lsr_num)
        cf: CF
        sg_type: 'carney', 'binomial'

        """
        self.set_freq(cf)

        self._hsr_num = anf_num[0]
        self._msr_num = anf_num[1]
        self._lsr_num = anf_num[2]


        # Outer/middle ear filter
        # Stapes velocity [Pa -> m/s]
        # BM module is in run()
        # IHC receptor potential


        if self._hsr_num != 0:
            self.ihc_hsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr.read_pars(par_dir("ihc_hsr_Sumner2002.par"))

            self.anf_hsr = self._generate_anf(sg_type, 1)
            dsam.connect(self.ihc_hsr, self.anf_hsr)

        if self._msr_num != 0:
            self.ihc_msr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr.read_pars(par_dir("ihc_msr_Sumner2002.par"))

            self.anf_msr = self._generate_anf(sg_type, 1)
            dsam.connect(self.ihc_msr, self.anf_msr)


        if self._lsr_num != 0:
            self.ihc_lsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr.read_pars(par_dir("ihc_lsr_Sumner2002.par"))

            self.anf_lsr = self._generate_anf(sg_type, 1)
            dsam.connect(self.ihc_lsr, self.anf_lsr)


    def set_freq(self, cf):

        if isinstance(cf, int):
            cf = float(cf)

        # Only real numbers please.
        assert (isinstance(cf, float) or
                cf == None)

        if isinstance(cf, float):
            real_freq_map = tw.bm_pars.real_freq_map
            assert cf in real_freq_map
            self._freq_idx = int(np.where(real_freq_map == cf)[0])
        elif cf == None:
            self._freq_idx = None


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
        xBM = tw.run_bm(fs, sound, mode='x')

        ### IHCRP
        ihcrp = tw.run_ihcrp(fs, xBM)
        if self._freq_idx != None:
            ihcrp = ihcrp[:,self._freq_idx]


        trains = []
        if self._hsr_num > 0:
            self.ihc_hsr.run(fs, ihcrp)
            tr = self._run_anf('hsr', self.anf_hsr, fs, self._hsr_num)
            trains.extend(tr)

        if self._msr_num > 0:
            self.ihc_msr.run(fs, ihcrp)
            tr = self._run_anf('msr', self.anf_msr, fs, self._msr_num)
            trains.extend(tr)

        if self._lsr_num > 0:
            self.ihc_lsr.run(fs, ihcrp)
            tr = self._run_anf('lsr', self.anf_lsr, fs, self._lsr_num)
            trains.extend(tr)

        trains = np.rec.array(trains, dtype=self._train_type)

        return trains


def main():
    import thorns as th

    fs = 48000
    cf = tw.real_freq_map[38]
    print "CF:", cf
    stimdb = 70

    ear = Holmberg2008((250,0,0), cf=cf)

    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = dsam.set_dbspl(stimdb, s)
    z = np.zeros( np.ceil(len(t)/4) )
    s = np.concatenate( (z, s, z) )


    anf = ear.run(fs, s)
    th.plot_raster(anf['spikes']).show()
    th.plot_psth(anf['spikes'], bin_size=1).show()


if __name__ == "__main__":
    main()
