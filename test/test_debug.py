import numpy as np
import matplotlib.pyplot as plt

import bai_bm

def main():

    print "Start:"

    debug_orig = np.fromfile("uihc.dat", float)
    debug_orig = debug_orig.reshape((len(debug_orig)/100, 100))
    debug_orig = np.fliplr(debug_orig)



    # Generate debug.dat
    fs = 48000                  # Hz
    xBM_bin = np.fromfile("xbm.dat", float)
    xBM_bin = xBM_bin.reshape((len(xBM_bin)/100, 100))
    # xBM_bin = np.fliplr(xBM_bin)
    uIHC = bai_bm.run_ihcrp(fs, xBM_bin)


    # Read debug.dat
    debug = np.fromfile("debug.dat", float)
    debug = debug.reshape((len(debug)/100, 100))
    debug = np.fliplr(debug)
    debug = np.fliplr(uIHC)

    # plt.imshow(xBM_target, aspect='auto')
    # plt.colorbar()
    # plt.show()

    # plt.imshow(xBM, aspect='auto')
    # plt.colorbar()
    # plt.show()

    # sec = 30
    # plt.plot(xBM_target[:,sec])
    # plt.plot(xBM[:,sec])
    # plt.show()

    print debug_orig.shape, debug.shape

    plt.imshow(debug_orig, aspect='auto')
    plt.show()
    plt.imshow(debug, aspect='auto')
    plt.show()


    sec = 60
    plt.plot(debug_orig[:,sec])
    plt.plot(debug[:,sec])
    plt.show()

    # plt.plot(debug[:,sec]-uIHC[:,sec])


    plt.imshow(np.abs(debug_orig[0:-1,:] - debug), aspect='auto')
    plt.colorbar()
    plt.show()


    print "done."

if __name__ == "__main__":
    main()
