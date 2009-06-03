import numpy as np
import scipy.signal as dsp
import matplotlib.pyplot as plt


def plot_mid_ear_filter():
    """
    Plot middle ear filter's characteristic from Werner's model.
    """
    fs = 48000.0
    S_ST = 3.1e-6                          # Area of stapes [m^2]
    S_ED = 55.e-6;                         # Ear drum area in m^2
    C_eardrum = (0.7e-9/20.e-3/S_ED);      # Ear drum compliance nm/dbspl/m^2


    # Digital wave filter coefficients
    R2_ME = 1. / (2. * fs * C_eardrum)
    R1 = 1. / (2. * np.pi * C_eardrum * 1e3)
    g1_ME = (R2_ME - R1) / (R2_ME + R1)
    Z_ME=0

    # Standard filter coefficients
    b = [(-1-g1_ME), (1+g1_ME)]
    a = [R2_ME, R2_ME*g1_ME]


    w, h = dsp.freqz(b, a)

    f = np.linspace(0, fs/2.0, len(w))


    plt.semilogx(f, 20*np.log10(abs(h)))
    plt.show()


if __name__ == "__main__":
    plot_mid_ear_filter()
