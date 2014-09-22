# -*- coding: utf-8 -*-

"""Interoperability with Brian.

"""

from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals

__author__ = "Marek Rudnicki"
__copyright__ = "Copyright 2014, Marek Rudnicki"
__license__ = "GPLv3+"


import numpy as np

import cochlea
import thorns as th
import thorns.waves as wv

def main():

    fs = 100e3
    cf = 1e3

    sound = wv.ramped_tone(
        fs=fs,
        freq=cf,
        duration=50e-3,
        pad=30e-3,
        dbspl=50,
    )

    anf_trains = cochlea.run_zilany2014(
        sound=sound,
        fs=fs,
        cf=cf,
        anf_num=(300, 0, 0),
        seed=0,
        species='cat',
    )

    anfs = cochlea.make_brian_group(anf_trains)

    print(anfs)


if __name__ == "__main__":
    main()
