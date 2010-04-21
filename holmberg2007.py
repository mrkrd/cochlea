from __future__ import division

import numpy as np

from auditory_periphery import par_dir, AuditoryPeriphery
import traveling_waves as tw
import dsam


class Holmberg2007(AuditoryPeriphery):
    def __init__(self, anf_num=(1,1,1), cf=None,
                 sg_type='carney', accumulate=False):
        """ Auditory periphery model from Marcus Holmberg

        anf_num: (hsr_num, msr_num, lsr_num)
        cf: CF
        sg_type: 'carney', 'binomial'
        accumulate: if True, spikes for all fibers are calculated at once

        """
        self.name = "Holmberg2007"

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

            self.sg_hsr_module = self._generate_anf(self._hsr_num, sg_type, accumulate)
            dsam.connect(self.ihc_hsr_module, self.sg_hsr_module)

        if self._msr_num > 0:
            self.ihc_msr_module = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr_module.read_pars(par_dir("ihc_msr_Sumner2002.par"))

            self.sg_msr_module = self._generate_anf(self._msr_num, sg_type, accumulate)
            dsam.connect(self.ihc_msr_module, self.sg_msr_module)


        if self._lsr_num > 0:
            self.ihc_lsr_module = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr_module.read_pars(par_dir("ihc_lsr_Sumner2002.par"))

            self.sg_lsr_module = self._generate_anf(self._lsr_num, sg_type, accumulate)
            dsam.connect(self.ihc_lsr_module, self.sg_lsr_module)


    def set_freq(self, cf):

        if isinstance(cf, int):
            cf = float(cf)

        # Only real numbers please.
        assert isinstance(cf, float) | (cf is None)

        if isinstance(cf, float):
            real_freq_map = tw.bm_pars.real_freq_map
            assert cf in real_freq_map
            self._freq_idx = int(np.where(real_freq_map == cf)[0])
        elif cf is None:
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
        if self._freq_idx is not None:
            ihcrp = ihcrp[:,self._freq_idx]


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

        trains = np.array(trains, dtype=self._train_type)

        return trains


def main():
    import thorns as th

    fs = 48000
    cf = tw.real_freq_map[38]
    print "CF:", cf
    stimdb = 70

    ear = Holmberg2007((250,0,0), cf=cf)

    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = dsam.set_dbspl(stimdb, s)
    z = np.zeros( np.ceil(len(t)/4) )
    s = np.concatenate( (z, s, z) )


    anf = ear.run(fs, s)
    th.plot_raster(anf['spikes']).show()
    th.plot_psth(anf['spikes'], bin_size=1).show()


    # ear = Holmberg2007((1,0,0))
    # anf = ear.run(fs, s)
    # th.plot_raster(anf['spikes']).show()
    # th.plot_psth(anf['spikes'], bin_size=1).show()



if __name__ == "__main__":
    main()
