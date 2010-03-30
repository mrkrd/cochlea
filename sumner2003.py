# Author: Marek Rudnicki
# Time-stamp: <2010-03-16 21:04:08 marek>
#
# Description: Sumner et al. ``A nonlinear filter-bank model of the
# guinea-pig cochlear nerve: Rate responses''
# J. Acoust. Soc. Am. Volume 113, Issue 6, pp. 3264-3274 (June 2003)

from __future__ import division

import numpy as np

import thorns as th
import dsam

from auditory_periphery import AuditoryPeriphery, par_dir

class Sumner2003(AuditoryPeriphery):
    def __init__(self, anf_num=(1, 1, 1), cf=1000,
                 sg_type='carney', accumulate=False):
        """ Auditory periphery model from Sumner (2003)

        anf_num: (hsr_num, msr_num, lsr_num)
        cf: CF
        sg_type: 'carney', 'binomial'
        accumulate: if True, then spike trains of each type are concatenated

        """
        self.name = 'Sumner2003'

        assert accumulate == False

        self._hsr_num = anf_num[0]
        self._msr_num = anf_num[1]
        self._lsr_num = anf_num[2]


        # TODO: check the filters
        self.outer_middle_ear_A = dsam.EarModule("Filt_MultiBPass")
        self.outer_middle_ear_A.read_pars(par_dir("filt_GPa_Sumner2003.par"))

        self.outer_middle_ear_B = dsam.EarModule("Filt_MultiBPass")
        self.outer_middle_ear_B.read_pars(par_dir("filt_GPb_Sumner2003.par"))
        dsam.connect(self.outer_middle_ear_A, self.outer_middle_ear_B)


        # Stapes velocity [Pa -> m/s]
        self.stapes_velocity = dsam.EarModule("Util_mathOp")
        self.stapes_velocity.read_pars(par_dir("stapes_Sumner2003.par"))
        dsam.connect(self.outer_middle_ear_B, self.stapes_velocity)


        # Basilar membrane
        self.bm = dsam.EarModule("BM_DRNL")
        self.bm.read_pars(par_dir("bm_Sumner2003.par"))
        self.set_freq(cf)
        dsam.connect(self.stapes_velocity, self.bm)


        # IHC receptor potential
        self.ihcrp = dsam.EarModule("IHCRP_Shamma3StateVelIn")
        self.ihcrp.read_pars(par_dir("ihcrp_Sumner2002.par"))
        dsam.connect(self.bm, self.ihcrp)

        if self._hsr_num > 0:
            self.ihc_hsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr.read_pars(par_dir("ihc_hsr_Sumner2003.par"))
            dsam.connect(self.ihcrp, self.ihc_hsr)

            self.anf_hsr = self._generate_anf(sg_type, 1)
            dsam.connect(self.ihc_hsr, self.anf_hsr)

        if self._msr_num > 0:
            self.ihc_msr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr.read_pars(par_dir("ihc_msr_Sumner2003.par"))
            dsam.connect(self.ihcrp, self.ihc_msr)

            self.anf_msr = self._generate_anf(sg_type, 1)
            dsam.connect(self.ihc_msr, self.anf_msr)


        if self._lsr_num > 0:
            self.ihc_lsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr.read_pars(par_dir("ihc_lsr_Sumner2003.par"))
            dsam.connect(self.ihcrp, self.ihc_lsr)

            self.anf_lsr = self._generate_anf(sg_type, 1)
            dsam.connect(self.ihc_lsr, self.anf_lsr)



    def set_freq(self, cf):
        # Test for `cf' type, should be either float or tuple
        # flaot: single BM frequency
        # tuple: (min_freq, max_freq, freq_num)
        if isinstance(cf, int):
            cf = float(cf)
        assert (isinstance(cf, tuple) or
                isinstance(cf, float))

        if isinstance(cf, float):
            self.bm.set_par("CF_MODE", "single")
            self.bm.set_par("SINGLE_CF", cf)
        elif isinstance(cf, tuple):
            self.bm.set_par("CF_MODE", "guinea_pig")
            self.bm.set_par("MIN_CF", cf[0])
            self.bm.set_par("MAX_CF", cf[1])
            self.bm.set_par("CHANNELS", cf[2])


    def get_freq_map(self):
        """ Returns the current frequancy map.

        Note: since frequency map is created during simulation, this
        function hast to be called after run()

        """
        return self.bm.get_labels()


    def print_pars(self):
        self.bm.print_pars()
        if self._hsr_num > 0:
            self.ihc_hsr.print_pars()

        if self._msr_num > 0:
            self.ihc_msr.print_pars()

        if self._lsr_num > 0:
            self.ihc_lsr.print_pars()


    def run(self, fs, sound):
        """ Run auditory periphery model.

        fs: sampling frequency
        sound: audio signal

        """
        self.outer_middle_ear_A.run(fs, sound)
        self.outer_middle_ear_B.run()

        self.stapes_velocity.run()
        self.bm.run()
        self.ihcrp.run()

        trains = []
        if self._hsr_num > 0:
            self.ihc_hsr.run()
            tr = self._run_anf('hsr', self.anf_hsr, fs, self._hsr_num)
            trains.extend(tr)

        if self._msr_num > 0:
            self.ihc_msr.run()
            tr = self._run_anf('msr', self.anf_msr, fs, self._msr_num)
            trains.extend(tr)

        if self._lsr_num > 0:
            self.ihc_lsr.run()
            tr = self._run_anf('lsr', self.anf_lsr, fs, self._lsr_num)
            trains.extend(tr)

        trains = np.rec.array(trains, dtype=self._train_type)

        return trains



def main():
    fs = 100000
    cf = 1000
    stimdb = 70

    ear = Sumner2003((250,0,0), cf=cf)

    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = dsam.set_dbspl(stimdb, s)
    z = np.zeros(np.ceil(len(t)/3))
    s = np.concatenate( (z, s, z) )

    anf = ear.run(fs, s)

    p = th.plot_raster(anf['spikes'])
    p.show()

    p = th.plot_psth(anf['spikes'])
    p.show()


if __name__ == "__main__":
    main()

