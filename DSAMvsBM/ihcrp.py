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
    tone.set_par("DT", 1.0/fs)

    # DRNL basilar membrane
    stapes = dsam.EarModule("Util_mathOp")
    stapes.read_pars("stapes.par")
    stapes.print_pars()

    bm_drnl = dsam.EarModule("BM_DRNL")
    bm_drnl.read_pars("bm_drnl_human.par")
    bm_drnl.print_pars()

    ihcrp = dsam.EarModule("IHCRP_Shamma3StateVelIn")
    ihcrp.read_pars("ihcrp_Meddis2005.par")
    ihcrp.print_pars()


    dsam.connect(tone, stapes)
    dsam.connect(stapes, bm_drnl)
    dsam.connect(bm_drnl, ihcrp)

    tone.run()
    stapes.run()
    bm_drnl.run()
    ihcrp.run()

    ihcrp_dsam = ihcrp.get_signal()

    plt.imshow(ihcrp_dsam.T, aspect='auto')
    plt.colorbar()
    plt.show()


    # BAI basilar membrane
    tone_arr = tone.get_signal()
    tone_arr = bai_bm.run_stapes(tone_arr)
    bm_bai_x = bai_bm.run_bm(fs, tone_arr, mode='x')
    ihcrp_bai = bai_bm.run_ihcrp(fs, bm_bai_x)

    plt.imshow(ihcrp_bai.T, aspect='auto')
    plt.colorbar()
    plt.show()

    print "Max indexes:"
    print "DSAM: ", np.sum(ihcrp_dsam, axis=0).argmax()
    print " BAI: ", np.sum(ihcrp_bai, axis=0).argmax()

    plt.plot(ihcrp_dsam[:,38])
    plt.plot(ihcrp_bai[:,38])
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
