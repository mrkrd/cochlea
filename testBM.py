import sys
import numpy as np
import matplotlib.pyplot as plt
from bai_bm import bm



signal = np.load('test/forward.npy')


def main():


    xBM = bm.run(signal)

    orig_xBM = np.load('test/xBM.npy')
    orig_LCR4 = np.load('test/LCR4.npy')

    d = xBM - orig_LCR4
    print np.max(xBM), np.max(orig_LCR4), np.max(d)

    plt.plot(xBM[:,75])
    plt.plot(orig_LCR4[:,75])
    plt.show()



    # plt.imshow(xBM.T, aspect='auto')

    # for i in range(5):
    #     plt.plot(xBM[:,i*20])

    #plt.plot(xBM[:,70])
    # plt.imshow(xBM.T, aspect='auto')
    # plt.show()


if __name__ == "__main__":
    import cProfile
    # for i in range(1000):
    #     print i
    #     main()

    main()

    # cProfile.run('main()')
