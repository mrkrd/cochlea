from __future__ import division

import numpy as np

from zilany2009 import (
    Zilany2009,
    run_zilany2009
)

from zilany2009_human import (
    Zilany2009_Human,
    run_zilany2009_human
)

from zilany2009_human_psp import Zilany2009_Human_PSP

try:
    from holmberg2007 import Holmberg2007
except ImportError:
    print "Holmberg2007 not loaded (DSAM missing?)"




def set_dbspl(signal, dbspl):
    p0 = 20e-6
    rms = np.sqrt( np.sum(signal**2) / signal.size )

    scalled = signal * 10**(dbspl/20) * p0 / rms

    return scalled


def set_dba_isolet(signal, dba):
    p0 = 20e-6
    rms_dba = 0.02972401089     # value from miclib

    scaled = signal * 10**(dba/20) * p0 / rms_dba

    return scaled
