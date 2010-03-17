import numpy as np
import biggles
from scipy.io import loadmat
import sys

import traveling_waves as tw


def plot_freq_maps(files):
    """
    Plot real_freq_map of BM model against all freq_maps found in
    matlab files.
    """

    lin = loadmat('../bm_pars_orig/param_lin.mat', squeeze_me=True)
    res = loadmat('../bm_pars_orig/param_res.mat', squeeze_me=True)

    freq_map_greenwood = lin['freq_map_greenwood'][::-1]
    freq_map = lin['freq_map'][::-1]
    freq_map_wished = lin['freq_map_wished'][::-1]
    freq_map_res = res['freq_map_res'][::-1]

    # real_freq_map_file = np.load('real_freq_map.npy')


    secs = np.arange(100)

    p = biggles.FramedPlot()
    # p.ylog = 1
    p.add( biggles.Curve(secs, freq_map_greenwood, type='dashed', width=3) )
    p.add( biggles.Curve(secs, freq_map) )
    p.add( biggles.Curve(secs, freq_map_wished) )
    p.add( biggles.Curve(secs, freq_map_res) )

    for f in files:
        m = np.load(f)
        p.add( biggles.Curve(secs, m, color='red') )

    p.write_eps("freq_maps.eps")
    p.show()


if __name__ == "__main__":
    files = sys.argv[1:]
    plot_freq_maps(files)
