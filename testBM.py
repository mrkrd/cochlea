import sys
import numpy as np
import matplotlib.pyplot as plt
import _bm




bm_pars = np.load('bm_pars.npz')

fs = 48000.0
t = np.arange(0, 0.1, 1/fs)
signal = np.zeros(len(t))
signal[len(signal)/5] = 1e-9;


signal = np.load('test/forward.npy')


def main():

    _bm.bm_init(48000,
                bm_pars['Ls'],
                bm_pars['Rs'],
                bm_pars['Ct'],
                bm_pars['Rbm'],
                bm_pars['Cbm'],
                bm_pars['Lbm'],
                float(bm_pars['Rh']),
                float(bm_pars['Lh']))


    xBM = _bm.bm_wave(signal,
                      bm_pars['ampl_corr'],
                      bm_pars['Abm'],
                      bm_pars['Cbm'])



    orig_xBM = np.load('test/xBM.npy')

    d = xBM - orig_xBM
    print np.max(xBM), np.max(orig_xBM), np.max(d)

    # plt.plot(xBM[:,75])
    # plt.plot(orig_xBM[:,75])
    # plt.show()


    # print sys.getrefcount(xBM)

    _bm.LCR4_init(fs,
                  bm_pars['freq_map'],
                  bm_pars['Qmin'],
                  bm_pars['SAT1'],
                  bm_pars['SAT4'])


    xBM = _bm.LCR4(xBM,
                   bm_pars['Qmax'],
                   bm_pars['Qmin']);

    # print sys.getrefcount(xBM)


    orig_LCR4 = np.load('test/LCR4.npy')

    d = xBM - orig_LCR4
    print np.max(xBM), np.max(orig_LCR4), np.max(d)

    plt.plot(xBM[:,75])
    plt.plot(orig_LCR4[:,75])
    plt.show()

    vBM = np.diff(xBM, axis=0) * fs


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