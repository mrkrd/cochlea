#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np

import dsam

def main():
    ihc = dsam.EarModule("IHC_Meddis2000")
    ihc.read_pars("ihc_hsr_Sumner2002.par")

    fs = 48e3

    s = np.loadtxt(
        "ihcrp_fs_48e3.out",
    )

    ihc.run(fs, s)
    r = ihc.get_signal()

    # import marlib as mr
    # mr.plot(r)

    np.savetxt(
        "ihc_meddis2000_dsam.out",
        r,
        fmt='%.70e'
    )


if __name__ == "__main__":
    main()
