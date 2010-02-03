# Author: Marek Rudnicki
# Time-stamp: <2010-01-30 18:29:08 marek>
#
# Description: Model of auditory periphery as described by Sumner et
# al. (2002)

from __future__ import division

import numpy as np

import thorns as th
import dsam

from auditory_periphery import AuditoryPeriphery, par_dir

class Sumner2002(AuditoryPeriphery):
    def __init__(self, anf_num=(1, 1, 1), freq=1000,
                 animal='gp', sg_type='carney', accumulate=False):
        """ Auditory periphery model from Sumner (2002)

        anf_num: (hsr_num, msr_num, lsr_num)
        freq: CF
        animal: 'gp', 'human'
        powerlaw_implnt: 'approx' or 'acctual' implementation of the power-law
        accumulate: if True, then spike trains of each type are concatenated

        """
        assert accumulate == False

        self._hsr_num = anf_num[0]
        self._msr_num = anf_num[1]
        self._lsr_num = anf_num[2]
        self._animal = animal


        # TODO: check the filters
        if self._animal == 'gp':
            self.outer_middle_ear = dsam.EarModule("Filt_MultiBPass")
            self.outer_middle_ear.read_pars(par_dir("filt_GP_A.par"))

            self.outer_middle_ear_B = dsam.EarModule("Filt_MultiBPass")
            self.outer_middle_ear_B.read_pars(par_dir("filt_GP_B.par"))
            dsam.connect(self.outer_middle_ear, self.outer_middle_ear_B)
        elif self._animal == 'human':
            self.outer_middle_ear = dsam.EarModule("Filt_MultiBPass")
            self.outer_middle_ear.read_pars(par_dir("filt_Human.par"))
        else:
            assert False


        # Stapes velocity [Pa -> m/s]
        self.stapes_velocity = dsam.EarModule("Util_mathOp")
        if self._animal == 'gp':
            self.stapes_velocity.read_pars(par_dir("stapes_Sumner2002.par"))
            dsam.connect(self.outer_middle_ear_B, self.stapes_velocity)
        elif self._animal == 'human':
            # TODO: fix the threshold for human
            self.stapes_velocity.set_par("OPERATOR", "SCALE")
            self.stapes_velocity.set_par("OPERAND", 1.7e-11)
            dsam.connect(self.outer_middle_ear, self.stapes_velocity)
        else:
            assert False


        # Basilar membrane
        if self._animal == 'gp':
            self.bm = dsam.EarModule("BM_DRNL")
            # TODO: fix BM parameters in order to have correct
            # rate-intensity curves (Nuria's work)
            self.bm.read_pars(par_dir("bm_drnl_gp.par"))
            self.set_freq(freq)
            dsam.connect(self.stapes_velocity, self.bm)
        elif self._animal == 'human':
            self.bm = dsam.EarModule("BM_DRNL")
            self.bm.read_pars(par_dir("drnl_human_Lopez-Poveda2001.par"))
            self.set_freq(freq)
            dsam.connect(self.stapes_velocity, self.bm)
        else:
            assert False


        # IHC receptor potential
        self.ihcrp = dsam.EarModule("IHCRP_Shamma3StateVelIn")
        self.ihcrp.read_pars(par_dir("ihcrp_Meddis2005_modified.par"))
        dsam.connect(self.bm, self.ihcrp)

        if self._hsr_num > 0:
            self.ihc_hsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr.read_pars(par_dir("ihc_hsr_Meddis2002.par"))
            dsam.connect(self.ihcrp, self.ihc_hsr)

            self.anf_hsr = self._generate_anf(sg_type, 1)
            dsam.connect(self.ihc_hsr, self.anf_hsr)

        if self._msr_num > 0:
            self.ihc_msr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr.read_pars(par_dir("ihc_msr_Meddis2002.par"))
            dsam.connect(self.ihcrp, self.ihc_msr)

            self.anf_msr = self._generate_anf(sg_type, 1)
            dsam.connect(self.ihc_msr, self.anf_msr)


        if self._lsr_num > 0:
            self.ihc_lsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr.read_pars(par_dir("ihc_lsr_Meddis2002.par"))
            dsam.connect(self.ihcrp, self.ihc_lsr)

            self.anf_lsr = self._generate_anf(sg_type, 1)
            dsam.connect(self.ihc_lsr, self.anf_lsr)



    def set_freq(self, freq):
        # Test for `freq' type, should be either float or tuple
        # flaot: single BM frequency
        # tuple: (min_freq, max_freq, freq_num)
        if isinstance(freq, int):
            freq = float(freq)
        assert (isinstance(freq, tuple) or
                isinstance(freq, float))

        if isinstance(freq, float):
            self.bm.set_par("CF_MODE", "single")
            self.bm.set_par("SINGLE_CF", freq)
        elif isinstance(freq, tuple):
            if self._animal == 'gp':
                self.bm.set_par("CF_MODE", "guinea_pig")
            elif self._animal == 'human':
                self.bm.set_par("CF_MODE", "human")
            else:
                assert False
            self.bm.set_par("MIN_CF", freq[0])
            self.bm.set_par("MAX_CF", freq[1])
            self.bm.set_par("CHANNELS", freq[2])


    def get_freq_map(self):
        """ Returns the current frequancy map.

        Note: since frequency map is created during simulation, this
        function hast to be called after run()

        """
        return self.bm.get_labels()


    def run(self, fs, sound):
        """ Run auditory periphery model.

        fs: sampling frequency
        sound: audio signal

        """
        self.outer_middle_ear.run(fs, sound)
        if self._animal == 'gp':
            self.outer_middle_ear_B.run()

        self.stapes_velocity.run()
        self.bm.run()
        self.ihcrp.run()

        hsr_trains = []
        msr_trains = []
        lsr_trains = []

        if self._hsr_num > 0:
            self.ihc_hsr.run()
            hsr_trains = self._run_anf(self.anf_hsr, fs, self._hsr_num)

        if self._msr_num > 0:
            self.ihc_msr.run()
            msr_trains = self._run_anf(self.anf_msr, fs, self._msr_num)

        if self._lsr_num > 0:
            self.ihc_lsr.run()
            lsr_trains = self._run_anf(self.anf_lsr, fs, self._lsr_num)

        hsr_trains = np.array(hsr_trains, dtype=self._train_type)
        msr_trains = np.array(msr_trains, dtype=self._train_type)
        lsr_trains = np.array(lsr_trains, dtype=self._train_type)

        return hsr_trains, msr_trains, lsr_trains



def main():
    ear = Sumner2002((250,0,0))

    fs = 100000
    cf = 1000
    stimdb = 50

    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * cf)
    s = dsam.set_dbspl(stimdb, s)
    z = np.zeros(np.ceil(len(t)/2))
    s = np.concatenate( (z, s, z) )

    ear.set_freq( cf )

    hsr, msr, lsr = ear.run(fs, s)
    th.plot_raster(hsr['spikes'])
    th.plot_psth(hsr['spikes'])


if __name__ == "__main__":
    main()

