import numpy as np
import matplotlib.pyplot as plt

import dsam
import bai_bm

def main():

    fs = 48000                  # Hz

    forward = np.load("forward.npy")
    forward = forward * 1e6     # Pa -> uPa

    ihcrp_target = np.load("ihcrp.npy")
    ihcrp_target = np.fliplr(ihcrp_target)

    bm_v = bai_bm.run_bm(fs, forward, mode='v')
    # bm_v = bm_v / 2.8

    dsam_input = dsam.EarModule(fs, bm_v)
    ihcrp = dsam.EarModule("IHCRP_Shamma3StateVelIn")
    ihcrp.read_pars("ihcrp_Meddis2005.par")
    dsam.connect(dsam_input, ihcrp)
    ihcrp.run()
    ihcrp_u = ihcrp.get_signal()


    sec = 61
    plt.plot(ihcrp_target[:,sec])
    plt.plot(ihcrp_u[:,sec])
    # plt.imshow(ihcrp_target, aspect='auto')
    # plt.colorbar()
    # plt.show()

    # plt.imshow(ihcrp_u, aspect='auto')
    # plt.colorbar()
    plt.show()



if __name__ == "__main__":
    main()
