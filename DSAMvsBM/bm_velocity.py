import sys
import numpy as np
import matplotlib.pyplot as plt
import bai_bm
import dsam



def main():
    fs = 48000.0

    tone = dsam.EarModule("Stim_PureTone_2")
    tone.set_par("FREQUENCY", 1000)
    tone.set_par("INTENSITY", 50)
    tone.set_par("DT", 1/fs)

    # DRNL basilar membrane
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

    plt.imshow(bm_dsam_v.T, aspect='auto')
    plt.colorbar()
    plt.show()


    # BAI basilar membrane
    tone_arr = tone.get_signal()
    tone_arr = tone_arr * bai_bm.S_ST * bai_bm.S_ED
    bm_bai_v = bai_bm.run_bm(fs, tone_arr, mode='v')

    plt.imshow(bm_bai_v.T, aspect='auto')
    plt.colorbar()
    plt.show()


    plt.plot(bm_dsam_v[:,61])
    plt.plot(bm_bai_v[:,61])
    plt.show()


    # plt.imshow(xBM.T, aspect='auto')

    # for i in range(5):
    #     plt.plot(xBM[:,i*20])

    # plt.plot(xBM[:,70])
    # plt.plot(vBM[:,70])
    # plt.imshow(np.rot90(vBM), aspect='auto')
    # plt.show()


if __name__ == "__main__":
    #import cProfile
    # for i in range(1000):
    #     print i
    #     main()

    main()

    # cProfile.run('main()')
