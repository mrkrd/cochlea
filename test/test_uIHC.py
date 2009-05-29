import numpy as np
import matplotlib.pyplot as plt

import bai_bm

def main():

    print "Start:"

    uIHC_target = np.fromfile("uihc.dat", float)
    uIHC_target = uIHC_target.reshape((len(uIHC_target)/100, 100))
    uIHC_target = np.fliplr(uIHC_target)


    fs = 48000                  # Hz
    xBM_bin = np.fromfile("xbm.dat", float)
    xBM_bin = xBM_bin.reshape((len(xBM_bin)/100, 100))
    xBM_bin = np.fliplr(xBM_bin)
    uIHC = bai_bm.run_ihcrp(fs, xBM_bin)


    plt.imshow(uIHC_target, aspect='auto')
    plt.colorbar()
    plt.show()

    plt.imshow(uIHC, aspect='auto')
    plt.colorbar()
    plt.show()

    plt.imshow(abs(uIHC_target-uIHC), aspect='auto')
    plt.colorbar()
    plt.show()

    # sec = 30
    # plt.plot(xBM_target[:,sec])
    # plt.plot(xBM[:,sec])
    # plt.show()

    sec = 30
    plt.plot(uIHC_target[:,sec])
    plt.plot(uIHC[:,sec])
    # plt.plot(uIHC_target[:,sec]-uIHC[:,sec])
    plt.show()

    print "done."

if __name__ == "__main__":
    main()
