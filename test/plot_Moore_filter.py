import numpy as np
import scipy.signal as dsp
import scipy.signal.filter_design as filt
import matplotlib.pyplot as plt
import sys

import dsam

# def dBfactor(dB, signal):
#     p0 = 2e-5                   # Pa
#     squared = signal ** 2
#     rms = np.sqrt( np.sum(squared) / len(signal) )

#     if (rms == 0):
#         r = 0.0
#     else:
#         r = 10**(dB / 20) * p0 / rms;

#     return r



def plot_Moore_filter():

    #b, a = filt.iirfilter(1, [0.3333333333, 0.666666666666])
    #w, h = dsp.freqz(b, a)

    # Input
    fs = 100000.0
    t = np.arange(0, 0.1, 1.0/fs)

    freq_list = np.arange(100, 20000, 100)

    stim_arr = np.zeros( (len(t), len(freq_list)) )

    for i,f in enumerate(freq_list):
        stim_arr[:,i] = np.sin(2 * np.pi * t * f)


    stim = dsam.EarModule(fs, stim_arr)

    # Filter
    filt = dsam.EarModule("Filt_MultiBPass")
    filt.read_pars("filt_moore.par")

    setdb = dsam.EarModule("Trans_SetDBSPL")
    setdb.set_par("dbspl", 0)



    intensity = dsam.EarModule("Ana_Intensity")


    # Connect DSAM objects
    dsam.connect(stim, setdb)
    dsam.connect(setdb, filt)
    dsam.connect(filt, intensity)

    setdb.run()
    filt.run()
    intensity.run()

    intensity_arr = intensity.get_signal()

    plt.semilogx(freq_list, intensity_arr)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplification (dB)")
    # plt.plot(intensity_arr)
    # plt.plot(filt.get_signal())
    plt.show()



if __name__ == "__main__":
    plot_Moore_filter()
