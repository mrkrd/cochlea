import numpy as np
import matplotlib.pyplot as plt

import dsam
import bai_bm

def main():

    fs = 48000                  # Hz

    forward = np.load("forward.npy")
    forward = forward * 1e6     # Pa -> uPa

    xBM_target = np.load("LCR4.npy")
    xBM_target = np.fliplr(xBM_target)

    bm_x = bai_bm.run_bm(fs, forward, mode='x')

    sec = 70
    plt.plot(xBM_target[:,sec])
    plt.plot(bm_x[:,sec])
    # plt.imshow(ihcrp_target, aspect='auto')
    # plt.colorbar()
    # plt.show()

    # plt.imshow(ihcrp_u, aspect='auto')
    # plt.colorbar()
    plt.show()



if __name__ == "__main__":
    main()
