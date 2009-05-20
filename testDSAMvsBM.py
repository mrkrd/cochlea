import sys
import numpy as np
import matplotlib.pyplot as plt
import _bm

sys.path.append('../pyDSAM')
import dsam


bm_pars = np.load('bm_pars.npz')

fs = 48000.0
t = np.arange(0, 0.1, 1/fs)
signal = np.zeros(len(t))
signal[len(signal)/5] = 1e-9


def main():

    tone = dsam.EarModule("Stim_PureTone_2")
    tone.set_par("FREQUENCY", 1000)
    tone.set_par("INTENSITY", 50)
    tone.set_par("DT", 1/fs)
    tone.print_pars()

    stapes = dsam.EarModule("Util_mathOp")
    stapes.read_pars("stapes.par")
    stapes.print_pars()

    bm_drnl = dsam.EarModule("BM_DRNL")
    bm_drnl.read_pars("bm_drnl_human.par")
    bm_drnl.print_pars()

    dsam.connect(tone, stapes)
    dsam.connect(stapes, bm_drnl)

    tone.run()
    stapes.run()
    bm_drnl.run()

    bm_dsam_v = bm_drnl.get_signal()

    # plt.imshow(bm_dsam_v.T, aspect='auto')
    # plt.show()


    _bm.bm_init(48000,
                bm_pars['Ls'],
                bm_pars['Rs'],
                bm_pars['Ct'],
                bm_pars['Rbm'],
                bm_pars['Cbm'],
                bm_pars['Lbm'],
                float(bm_pars['Rh']),
                float(bm_pars['Lh']))


    xBM = _bm.bm_wave(stapes.get_signal(),
                      bm_pars['ampl_corr'],
                      bm_pars['Abm'],
                      bm_pars['Cbm'])

    # print sys.getrefcount(xBM)

    _bm.LCR4_init(fs,
                  bm_pars['freq_map'],
                  bm_pars['Qmin'],
                  bm_pars['SAT1'],
                  bm_pars['SAT4'])


    xBM = _bm.LCR4(xBM,
                   bm_pars['Qmin'],
                   bm_pars['Qmax']);

    # print sys.getrefcount(xBM)


    vBM = np.diff(xBM, axis=0) * fs


    # plt.imshow(xBM.T, aspect='auto')

    # for i in range(5):
    #     plt.plot(xBM[:,i*20])

    # plt.plot(xBM[:,70])
    # plt.plot(vBM[:,70])
    # plt.imshow(np.rot90(vBM), aspect='auto')
    # plt.show()


if __name__ == "__main__":
    import cProfile
    # for i in range(1000):
    #     print i
    #     main()

    main()

    # cProfile.run('main()')
