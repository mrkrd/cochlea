import numpy as np
import os

import thorns as th
import dsam

from auditory_periphery import par_dir, AuditoryPeriphery

import traveling_waves as tw



class Holmberg2008(AuditoryPeriphery):
    def __init__(self, hsr=100, msr=100, lsr=100, freq=None, animal='human'):

        assert animal == 'human'

        self.set_freq(freq)

        self.hsr = hsr
        self.msr = msr
        self.lsr = lsr
        self.animal = animal

        # Outer/middle ear filter
        self.outer_middle_ear = dsam.EarModule("Filt_MultiBPass")
        self.outer_middle_ear.read_pars(par_dir("filt_Human.par"))


        # Stapes velocity [Pa -> m/s]
        self.stapes_velocity = dsam.EarModule("Util_mathOp")
        self.stapes_velocity.set_par("OPERATOR", "SCALE")
        # For Holmberg IHCRP module: OPERAND =  tw.S_ED * tw.S_ST * 13.5
        self.stapes_velocity.set_par("OPERAND", tw.S_ST * tw.S_ED * 3.3)
        dsam.connect(self.outer_middle_ear, self.stapes_velocity)


        # BM module is in run()


        # IHC receptor potential
        self.ihcrp = dsam.EarModule("IHCRP_Shamma3StateVelIn")
        self.ihcrp.read_pars(par_dir("ihcrp_Meddis2005_modified.par"))


        if self.hsr != 0:
            self.ihc_hsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr.read_pars(par_dir("ihc_hsr_Meddis2002.par"))
            dsam.connect(self.ihcrp, self.ihc_hsr)

            self.anf_hsr = dsam.EarModule("An_SG_Carney")
            self.anf_hsr.read_pars(par_dir("anf_carney.par"))
            self.anf_hsr.set_par("NUM_FIBRES", hsr)
            dsam.connect(self.ihc_hsr, self.anf_hsr)

        if self.msr != 0:
            self.ihc_msr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr.read_pars(par_dir("ihc_msr_Meddis2002.par"))
            dsam.connect(self.ihcrp, self.ihc_msr)

            self.anf_msr = dsam.EarModule("An_SG_Carney")
            self.anf_msr.read_pars(par_dir("anf_carney.par"))
            self.anf_msr.set_par("NUM_FIBRES", msr)
            dsam.connect(self.ihc_msr, self.anf_msr)


        if self.lsr != 0:
            self.ihc_lsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr.read_pars(par_dir("ihc_lsr_Meddis2002.par"))
            dsam.connect(self.ihcrp, self.ihc_lsr)

            self.anf_lsr = dsam.EarModule("An_SG_Carney")
            self.anf_lsr.read_pars(par_dir("anf_carney.par"))
            self.anf_lsr.set_par("NUM_FIBRES", lsr)
            dsam.connect(self.ihc_lsr, self.anf_lsr)


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
        """
        Returns frequency map of the model.
        """
        if isinstance(self._freq_idx, int):
            freq_map = tw.bm_pars.real_freq_map[self._freq_idx]
        else:
            freq_map = tw.bm_pars.real_freq_map

        return freq_map


    def run(self, fs, sound, times=1, output_format='spikes'):
        """
        Run auditory periphery model.

        sound: audio signal
        fs: sampling frequency
        times: how many many trials
        output_format: format of the output 'spikes' (for spiking
        times), 'signals' (for time function)
        """
        if output_format == 'signals':
            assert times == 1

        fs = float(fs)
        input_module = dsam.EarModule(fs, sound)

        dsam.connect(input_module, self.outer_middle_ear)
        self.outer_middle_ear.run()
        dsam.disconnect(input_module, self.outer_middle_ear)

        self.stapes_velocity.run()

        ### Basilar membrane
        filtered_signal = self.stapes_velocity.get_signal()
        bm_signal = tw.run_bm(fs, filtered_signal, mode='v')
        if self._freq_idx != None:
            bm_signal = bm_signal[:,self._freq_idx]


        ### IHCRP
        self.ihcrp.run(fs, bm_signal)

        # ihcrp_signal = tw.run_ihcrp(fs, bm_signal)
        # if self._freq_idx != None:
        #     ihcrp_signal = ihcrp_signal[:,self._freq_idx]
        # ihcrp_mod = dsam.EarModule(fs, ihcrp_signal)


        if self.hsr > 0:
            self.ihc_hsr.run()
            hsr_db = self._run_anf(self.anf_hsr, fs, times, output_format)
        else:
            hsr_db = None

        if self.msr > 0:
            self.ihc_msr.run()
            msr_db = self._run_anf(self.anf_msr, fs, times, output_format)
        else:
            msr_db = None

        if self.lsr > 0:
            self.ihc_lsr.run()
            lsr_db = self._run_anf(self.anf_lsr, fs, times, output_format)
        else:
            lsr_db = None


        return hsr_db, msr_db, lsr_db

