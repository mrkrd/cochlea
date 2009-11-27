from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def main():
    dat = np.load('si.npz')

    freq_range = dat['freq_range']
    si_sumner = dat['si_sumner']
    si_carney = dat['si_carney']


    # plt.pcolor(si_sumner)
    # plt.show()
    # plt.pcolor(si_holmberg)
    # plt.show()
    # plt.pcolor(si_carney)
    # plt.show()

    # data = plt.imread('si.png')

    fig = plt.gcf()
    ax = fig.add_subplot(111)

    # bg = ax.twinx()
    # bg.imshow(data)

    # scores_max_sumner = np.max(si_sumner, axis=1)
    # ax.semilogx(freq_range, scores_max_sumner)
    # scores_max_holmberg = np.max(si_holmberg, axis=1)
    # plt.semilogx(freq_range, scores_max_holmberg)
    scores_max_carney = np.max(si_carney, axis=1)
    ax.semilogx(freq_range, scores_max_carney)


    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Vector Strength")

    fig.savefig('carney2009_si.pdf')

    plt.close()



if __name__ == "__main__":
    main()
