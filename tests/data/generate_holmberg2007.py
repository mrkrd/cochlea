#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np

import marlib as mr
import cochlea
import dsam

from cochlea.holmberg2007.traveling_waves import (
    run_outer_ear_filter,
    run_middle_ear_filter,
    scaling_factor,
    run_bm_wave,
    run_lcr4,
    run_ihcrp
)


def main():
    fs = 48e3
    sound = np.zeros(0.05*fs)
    sound[int(0.01*fs)] = 1


    ### Outer ear filter
    sound_oe = run_outer_ear_filter(sound, fs)

    ### Middle ear filter
    sound_me = run_middle_ear_filter(sound_oe, fs)

    ### Scaling
    sound_scaled = sound_me * scaling_factor

    ### Basilar membrane
    xbm = run_bm_wave(sound_scaled, fs)

    ### Amplification
    lcr4 = run_lcr4(xbm, fs)

    ### IHCRP
    ihcrp = run_ihcrp(lcr4, fs)


    ### IHC/Synapse DSAM module
    ihc_module = dsam.EarModule("IHC_Meddis2000")
    ihc_module.read_pars("ihc_hsr_Sumner2002.par")
    ihc_module.run(fs, ihcrp)
    psp = ihc_module.get_signal()




    channel = 20

    np.savez_compressed(
        "holmberg2007",
        fs=fs,
        channel=channel,
        sound=sound,
        sound_oe=sound_oe,
        sound_me=sound_me,
        sound_scaled=sound_scaled,
        xbm=xbm[:,channel],
        lcr4=lcr4[:,channel],
        ihcrp=ihcrp[:,channel],
        psp=psp[:,channel]
    )


    mr.plot(sound)
    mr.plot(sound_oe)
    mr.plot(sound_me)
    mr.plot(xbm[:,channel])
    mr.plot(lcr4[:,channel])
    mr.plot(ihcrp[:,channel])
    mr.plot(psp[:,channel])


if __name__ == "__main__":
    main()
