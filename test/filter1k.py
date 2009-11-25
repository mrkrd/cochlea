import numpy as np
import matplotlib.pyplot as plt

import traveling_waves as tw

def main():
    fs = 100000.0
    fstim = 1000.0
    t = np.arange(0, 0.1, 1.0/fs)
    s = np.sin(2 * np.pi * t * fstim)

    sout = tw.run_mid_ear_filter(fs, s)

    plt.plot(t, s)
    plt.plot(t, sout)
    plt.show()


if __name__ == "__main__":
    main()
