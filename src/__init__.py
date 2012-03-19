import numpy as np

from zilany2009 import Zilany2009
from zilany2009_vesicles import Zilany2009_Vesicles
from zilany2009_voltage_vesicles import Zilany2009_Voltage_Vesicles
from zilany2009_human import Zilany2009_Human
from zilany2009_human_psp import Zilany2009_Human_PSP

try:
    from holmberg2007 import Holmberg2007
except ImportError:
    print "Holmberg2007 not loaded (DSAM missing?)"


def set_dbspl(signal, db):
    p0 = 2e-5                   # Pa
    squared = signal**2
    rms = np.sqrt( np.sum(squared) / len(signal) )

    if rms == 0:
        r = 0
    else:
        r = 10**(db / 20.0) * p0 / rms

    return signal * r           # Pa


