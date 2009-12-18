# Author: Marek Rudnicki
# Time-stamp: <2009-12-18 20:11:18 marek>
#
# Description: Model of auditory periphery as described by Sumner et
# al. (2002)

import numpy as np
import os

import thorns as th
import dsam

from auditory_periphery import AuditoryPeriphery, par_dir

class Sumner2002(AuditoryPeriphery):
    def __init__(self, anf_sum=(100, 100, 100),
                 freq=1000, animal='gp', sg_type='carney'):
        """ Auditory periphery model from Sumner (2002)

        anf_sum: (hsr_sum, msr_sum, lsr_sum)
        freq: CF
        animal: 'gp', 'human'

        """
        self._hsr_sum = anf_sum[0]
        self._msr_sum = anf_sum[1]
        self._lsr_sum = anf_sum[2]
        self._animal = animal


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
            self.stapes_velocity.read_pars(par_dir("stapes_Meddis2005.par"))
            dsam.connect(self.outer_middle_ear_B, self.stapes_velocity)
        elif self._animal == 'human':
            self.stapes_velocity.set_par("OPERATOR", "SCALE")
            self.stapes_velocity.set_par("OPERAND", 1.7e-11)
            dsam.connect(self.outer_middle_ear, self.stapes_velocity)
        else:
            assert False


        # Basilar membrane
        if self._animal == 'gp':
            self.bm = dsam.EarModule("BM_DRNL")
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

        if self._hsr_sum != 0:
            self.ihc_hsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr.read_pars(par_dir("ihc_hsr_Meddis2002.par"))
            dsam.connect(self.ihcrp, self.ihc_hsr)

            self.anf_hsr = self._generate_anf(sg_type, self._hsr_sum)
            dsam.connect(self.ihc_hsr, self.anf_hsr)

        if self._msr_sum != 0:
            self.ihc_msr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr.read_pars(par_dir("ihc_msr_Meddis2002.par"))
            dsam.connect(self.ihcrp, self.ihc_msr)

            self.anf_msr = self._generate_anf(sg_type, self._msr_sum)
            dsam.connect(self.ihc_msr, self.anf_msr)


        if self._lsr_sum != 0:
            self.ihc_lsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr.read_pars(par_dir("ihc_lsr_Meddis2002.par"))
            dsam.connect(self.ihcrp, self.ihc_lsr)

            self.anf_lsr = self._generate_anf(sg_type, self._lsr_sum)
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
        """
        Returns curretn frequancy map.

        Note: since frequency map is created during simulation, this
        function can be called after run()
        """
        return self.bm.get_labels()


    def run(self, fs, sound, times=1):
        """ Run auditory periphery model.

        fs: sampling frequency
        sound: audio signal
        times: how many many trials

        """
        self.outer_middle_ear.run(fs, sound)
        if self._animal == 'gp':
            self.outer_middle_ear_B.run()

        self.stapes_velocity.run()
        self.bm.run()
        self.ihcrp.run()

        if self._hsr_sum > 0:
            self.ihc_hsr.run()
            hsr_trains = self._run_anf(self.anf_hsr, fs, times)
        else:
            hsr_trains = None

        if self._msr_sum > 0:
            self.ihc_msr.run()
            msr_trains = self._run_anf(self.anf_msr, fs, times)
        else:
            msr_trains = None

        if self._lsr_sum > 0:
            self.ihc_lsr.run()
            lsr_trains = self._run_anf(self.anf_lsr, fs, times)
        else:
            lsr_trains = None


        return hsr_trains, msr_trains, lsr_trains



def main():
    ear = Sumner2002((1,0,0), freq=5000)

    fs = 100000.0
    cf = 1000
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

