import sys
import numpy as np
import matplotlib.pyplot as plt

import bai_bm


signal = np.load('test/forward.npy')


def main():

    xBM = bai_bm.run_bm(48000, signal, mode='v')


    plt.imshow(xBM.T, aspect='auto')
    plt.colorbar()
    plt.show()
    exit()


    orig_xBM = np.load('test/xBM.npy')
    orig_LCR4 = np.load('test/LCR4.npy')


    d = xBM - orig_LCR4
    print np.max(xBM), np.max(orig_LCR4), np.max(d)

    plt.plot(xBM[:,75])
    plt.plot(orig_LCR4[:,75])
    # plt.imshow(xBM - orig_LCR4, aspect='auto')
    plt.show()



if __name__ == "__main__":
    import cProfile
    # for i in range(1000):
    #     print i
    #     main()

    main()

    # cProfile.run('main()')
