from __future__ import division

import numpy as np

from auditory_periphery import par_dir, AuditoryPeriphery
import traveling_waves as tw
import dsam



class Holmberg2008(AuditoryPeriphery):
    def __init__(self, hsr=100, msr=100, lsr=100,
                 freq=None, animal='human', sg_type='carney'):

        assert animal == 'human'

        self.set_freq(freq)

        self.hsr_num = hsr
        self.msr_num = msr
        self.lsr_num = lsr
        self.animal = animal

        # Outer/middle ear filter
        self.outer_middle_ear = dsam.EarModule("Filt_MultiBPass")
        self.outer_middle_ear.read_pars(par_dir("filt_Human.par"))


        # Stapes velocity [Pa -> m/s]
        # BM module is in run()
        # IHC receptor potential


        if self.hsr_num != 0:
            self.ihc_hsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr.read_pars(par_dir("ihc_hsr_Meddis2002.par"))

            self.anf_hsr = self._generate_anf(sg_type, self.hsr_num)
            dsam.connect(self.ihc_hsr, self.anf_hsr)

        if self.msr_num != 0:
            self.ihc_msr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr.read_pars(par_dir("ihc_msr_Meddis2002.par"))

            self.anf_msr = self._generate_anf(sg_type, self.msr_num)
            dsam.connect(self.ihc_msr, self.anf_msr)


        if self.lsr_num != 0:
            self.ihc_lsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr.read_pars(par_dir("ihc_lsr_Meddis2002.par"))

            self.anf_lsr = self._generate_anf(sg_type, self.lsr_num)
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


        if self.hsr_num > 0:
            self.ihc_hsr.run(fs, ihcrp)
            hsr_db = self._run_anf(self.anf_hsr, fs, times, output_format)
        else:
            hsr_db = None

        if self.msr_num > 0:
            self.ihc_msr.run(fs, ihcrp)
            msr_db = self._run_anf(self.anf_msr, fs, times, output_format)
        else:
            msr_db = None

        if self.lsr_num > 0:
            self.ihc_lsr.run(fs, ihcrp)
            lsr_db = self._run_anf(self.anf_lsr, fs, times, output_format)
        else:
            lsr_db = None


        return hsr_db, msr_db, lsr_db



if __name__ == "__main__":
    import matplotlib.pyplot as plt

    fs = 48000
    fstim = 1000
    t = np.arange(0, 0.1, 1/fs)
    s = np.sin(2 * np.pi * t * fstim)
    s = dsam.set_dB_SPL(60, s)
    s = s * np.hanning(len(s))

    # s = np.zeros_like( t )
    # s[np.ceil(len(s)/3)] = 1000

    ear = Holmberg2008(hsr=100, msr=0, lsr=0)
    hsr, msr, lsr = ear.run(fs, s, output_format='signals')

    fig = plt.gcf()
    ax = fig.add_subplot(211)
    ax.plot(t, s)

    ax = fig.add_subplot(212)
    ax.imshow(hsr.T, aspect='auto', interpolation=None)
    # ax.colorbar()
    plt.show()
