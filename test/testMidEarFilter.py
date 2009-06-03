import numpy as np
import matplotlib.pyplot as plt

import bai_bm

fs = 48000.0
t = np.arange(0, 0.1, 1/fs)
signal = np.zeros(len(t))
signal[len(signal)/5] = 1
signal = np.sin(2 * np.pi * t * 100)

out_orig = bai_bm.run_mid_ear_filter_orig(fs, signal)
out = bai_bm.run_mid_ear_filter(fs, signal)


plt.plot(out_orig)
plt.plot(out)
plt.show()
