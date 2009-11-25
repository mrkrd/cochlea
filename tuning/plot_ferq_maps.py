import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

import traveling_waves as tw


def plot_freq_maps():
    """
    Plot real_freq_map of BM model against all freq_maps found in mat
    files.
    """
    lin = loadmat('../bm_pars_orig/param_lin.mat', squeeze_me=True)
    res = loadmat('../bm_pars_orig/param_res.mat', squeeze_me=True)

    freq_map_greenwood = lin['freq_map_greenwood'][::-1]
    freq_map = lin['freq_map'][::-1]
    freq_map_wished = lin['freq_map_wished'][::-1]
    freq_map_res = res['freq_map_res'][::-1]

    real_freq_map_file = np.load('real_freq_map.npy')

    fig = plt.gcf()

    # ax1 = fig.add_subplot(211)
    # ax1.semilogy(tw.real_freq_map, label="Real freq map (0dBSPL)")
    # ax1.semilogy(real_freq_map_file, label="Real freq map (80dBSPL)")
    # ax1.semilogy(freq_map_greenwood, label="freq_map_greenwood")
    # ax1.semilogy(freq_map, label="freq_map")
    # ax1.semilogy(freq_map_wished, label="freq_map_wished")
    # ax1.semilogy(freq_map_res, label="freq_map_res")
    # ax1.set_ylabel("Frequency [Hz]")

    ax2 = fig.add_subplot(111)
    ax2.plot(tw.real_freq_map, label="Real freq map (0dBSPL)")
    ax2.plot(real_freq_map_file, label="Real freq map (80dBSPL)")
    ax2.plot(freq_map_greenwood, label="freq_map_greenwood")
    ax2.plot(freq_map, label="freq_map")
    ax2.plot(freq_map_wished, label="freq_map_wished")
    ax2.plot(freq_map_res, label="freq_map_res")
    ax2.set_ylabel("Frequency [Hz]")
    ax2.set_xlabel("Section number")
    # ax2.set_xlim(70, 100)
    leg = ax2.legend(loc='best')
    for t in leg.get_texts():
        t.set_fontsize('small')    # the legend text fontsize


    fig.savefig("freq_map.eps")
    plt.show()

if __name__ == "__main__":
    plot_freq_maps()
