import numpy as np
import matplotlib.pyplot as plt

from audiogram_jackson99 import human_audiogram
from traveling_waves import real_freq_map

def plot_th():
    th = np.load('th.npz')

    ax = plt.gca()

    ax.semilogx(real_freq_map, th['sumner_th'], label="Sumner")
    ax.semilogx(real_freq_map, th['holmberg_th'], label="Holmberg")

    # Get only relevant freq range from human_audiogram
    max_freq = real_freq_map.max() + 2000
    min_freq = real_freq_map.min()
    relevant_range = ((human_audiogram['freq'] > min_freq) &
                      (human_audiogram['freq'] < max_freq))

    ax.semilogx(human_audiogram['freq'][relevant_range],
                human_audiogram['threshold'][relevant_range],
                label="Jackson")

    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Intensity [dB SPL]")
    ax.legend()

    plt.savefig("th.eps")

    plt.show()



if __name__ == "__main__":
    plot_th()
