import numpy as np
import matplotlib.pyplot as plt

import dsam

class Sumner2002(object):
    # TODO: add psth=False
    def __init__(self, hsr=100, msr=100, lsr=100, freq=1000.0):

        # Test for `freq' type, should be either float or tuple
        # flaot: single BM frequency
        # tuple: (min_freq, max_freq, freq_num)
        if isinstance(freq, int):
            freq = float(freq)
        assert isinstance(freq, tuple) or isinstance(freq, float)

        self.hsr = hsr
        self.msr = msr
        self.lsr = lsr

        # Outer/middle ear filter
        self.outer_middle_ear = dsam.EarModule("Filt_MultiBPass")
        self.outer_middle_ear.read_pars("pars/filt_Moore.par")


        # Stapes velocity [Pa -> m/s]
        self.stapes_velocity = dsam.EarModule("Util_mathOp")
        self.stapes_velocity.read_pars("pars/stapes_Meddis2005.par")
        dsam.connect(self.outer_middle_ear, self.stapes_velocity)


        # Basilar membrane
        self.bm = dsam.EarModule("BM_DRNL")
        # TODO: change name to bm_human_MeddisXXXX.par
        self.bm.read_pars("pars/bm_drnl_human.par")
        if isinstance(freq, float):
            self.bm.set_par("CF_MODE", "single")
            self.bm.set_par("SINGLE_CF", freq)
        elif isinstance(freq, tuple):
            self.bm.set_par("CF_MODE", "human")
            self.bm.set_par("MIN_CF", freq[0])
            self.bm.set_par("MAX_CF", freq[1])
            self.bm.set_par("CHANNELS", freq[2])
        dsam.connect(self.stapes_velocity, self.bm)


        # IHC receptor potential
        self.ihcrp = dsam.EarModule("IHCRP_Shamma3StateVelIn")
        self.ihcrp.read_pars("pars/ihcrp_Meddis2005.par")
        dsam.connect(self.bm, self.ihcrp)

        if self.hsr != 0:
            self.ihc_hsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_hsr.read_pars("pars/ihc_hsr_Meddis2002.par")
            dsam.connect(self.ihcrp, self.ihc_hsr)

            self.anf_hsr = dsam.EarModule("An_SG_Binomial")
            self.anf_hsr.read_pars("pars/anf_binomial.par")
            self.anf_hsr.set_par("NUM_FIBRES", hsr)
            dsam.connect(self.ihc_hsr, self.anf_hsr)
            self.anf_hsr.print_pars()

        if self.msr != 0:
            self.ihc_msr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_msr.read_pars("pars/ihc_msr_Meddis2002.par")
            dsam.connect(self.ihcrp, self.ihc_msr)

            self.anf_msr = dsam.EarModule("An_SG_Binomial")
            self.anf_msr.read_pars("pars/anf_binomial.par")
            self.anf_msr.set_par("NUM_FIBRES", msr)
            dsam.connect(self.ihc_msr, self.anf_msr)


        if self.lsr != 0:
            self.ihc_lsr = dsam.EarModule("IHC_Meddis2000")
            self.ihc_lsr.read_pars("pars/ihc_lsr_Meddis2002.par")
            dsam.connect(self.ihcrp, self.ihc_lsr)

            self.anf_lsr = dsam.EarModule("An_SG_Binomial")
            self.anf_lsr.read_pars("pars/anf_binomial.par")
            self.anf_lsr.set_par("NUM_FIBRES", lsr)
            dsam.connect(self.ihc_lsr, self.anf_lsr)


    def run(self, fs, sound):
        fs = float(fs)

        input_module = dsam.EarModule(fs, sound)

        dsam.connect(input_module, self.outer_middle_ear)

        self.outer_middle_ear.run()
        self.stapes_velocity.run()
        self.bm.run()
        self.ihcrp.run()
        if self.hsr != 0:
            self.anf_hsr.set_par("PULSE_DURATION", 1/fs)
            self.ihc_hsr.run()
            self.anf_hsr.run()
            hsr_out = self.anf_hsr.get_signal()
        else:
            hsr_out = None

        if self.msr != 0:
            #self.anf_msr.set_par("PULSE_DURATION", 1/fs)
            self.ihc_msr.run()
            self.anf_msr.run()
            msr_out = self.anf_msr.get_signal()
        else:
            msr_out = None

        if self.lsr != 0:
            #self.anf_lsr.set_par("PULSE_DURATION", 1/fs)
            self.ihc_lsr.run()
            self.anf_lsr.run()
            lsr_out = self.anf_lsr.get_signal()
        else:
            lsr_out = None

        # plt.plot(self.ihcrp.get_signal())
        # plt.show()

        dsam.disconnect(input_module, self.outer_middle_ear)

        return hsr_out, msr_out, lsr_out





class Holmberg2008(object):
    def __init__(self, hsr=100, msr=100, lsr=100):

        self.outer_middle_ear = dsam.EarModule("Filt_MultiBPass")
        self.ihcrp = dsam.EarModule("IHCRP_Shamma3StateVelIn")
        self.ihc_hsr = dsam.EarModule("IHC_Meddis2000")
        self.ihc_msr = dsam.EarModule("IHC_Meddis2000")
        self.ihc_lsr = dsam.EarModule("IHC_Meddis2000")
        self.anf_hsr = dsam.EarModule("An_SG_Binomial")
        self.anf_msr = dsam.EarModule("An_SG_Binomial")
        self.anf_lsr = dsam.EarModule("An_SG_Binomial")

        self.outer_middle_ear.read_pars("pars/filt_Moore.par")
        self.ihcrp.read_pars("pars/ihcrp_Meddis2005.par")
        self.ihc_hsr.read_pars("pars/ihc_hsr_Meddis2002.par")
        self.ihc_msr.read_pars("pars/ihc_msr_Meddis2002.par")
        self.ihc_lsr.read_pars("pars/ihc_lsr_Meddis2002.par")
        self.anf_hsr.read_pars("pars/anf_binomial.par")
        self.anf_msr.read_pars("pars/anf_binomial.par")
        self.anf_lsr.read_pars("pars/anf_binomial.par")

        dsam.connect(self.ihcrp, self.ihc_hsr)
        dsam.connect(self.ihcrp, self.ihc_msr)
        dsam.connect(self.ihcrp, self.ihc_lsr)
        dsam.connect(self.ihc_hsr, self.anf_hsr)
        dsam.connect(self.ihc_msr, self.anf_msr)
        dsam.connect(self.ihc_lsr, self.anf_lsr)


    def run(self, fs, sound):
        self.anf_hsr.set_par("PULSE_DURATION", 0.000021)
        self.anf_msr.set_par("PULSE_DURATION", 0.000021)
        self.anf_lsr.set_par("PULSE_DURATION", 0.000021)


        # Filter input signal through outer/middle ear
        input_dsam = dsam.EarModule(fs, sound)
        dsam.connect(input_dsam, self.outer_middle_ear)
        self.outer_middle_ear.run()
        dsam.disconnect(input_dsam, self.outer_middle_ear)
        filtered = self.outer_middle_ear.get_signal()
        filtered = filtered * bai_bm.S_ED * bai_bm.S_ST # stapes velocity

        # Basilar membrane
        bm = bai_bm.run_bm(fs, filtered, mode='v')

        # IHCRP, IHC, SG
        bm_input_dsam = dsam.EarModule(fs, bm)
        dsam.connect(bm_input_dsam, self.ihcrp)
        self.ihcrp.run()
        dsam.disconnect(bm_input_dsam, self.ihcrp)
        self.ihc_hsr.run()
        self.ihc_msr.run()
        self.ihc_lsr.run()
        self.anf_hsr.run()
        self.anf_msr.run()
        self.anf_lsr.run()


        hsr_out = self.anf_hsr.get_signal()
        msr_out = self.anf_msr.get_signal()
        lsr_out = self.anf_lsr.get_signal()

        return hsr_out, msr_out, lsr_out






def main():
    import matplotlib.pyplot as plt

    print "start:"
    ear = Sumner2002(freq=(50,20000,100))

    fs = 100000.0
    t = np.arange(0, 0.1, 1.0/fs)
    s = np.sin(2 * np.pi * t * 1000)
    s = s * 3000                # ~40dBSPL
    s = s * np.hanning(len(s))


    hsr, msr, lsr = ear.run(fs, s)
    # hsr, msr, lsr = ear.run(fs, s)

    plt.imshow(hsr.T, aspect='auto', interpolation=None)
    plt.colorbar()
    plt.show()

    plt.plot(np.sum(hsr, axis=0))
    plt.show()

    m = np.sum(hsr, axis=0).argmax()
    plt.plot(hsr[:,m])
    plt.show()

    print "done."


if __name__ == "__main__":
    import cProfile

    # cProfile.run('main()')
    main()
