# Description: Model of auditory periphery of: Zilany, M.S.A., Bruce,
# I.C., Nelson, P.C., and Carney, L.H. (manuscript in preparation) 2009

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import scipy.signal as dsp
import os

from . import _pycat
from . zilany2009_human import (
    _calc_cfs,
    _run_human_me_filter_for_zilany2009
)


def run_zilany2009_human_psp(
        sound,
        fs,
        anf_type,
        cf,
        cohc=1,
        cihc=1,
        powerlaw_implnt='approx',
        with_ffGn=False):



    assert np.max(sound) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1


    cfs = _calc_cfs(cf)


    # Run Outer/Middle Ear filter
    sound = _run_human_me_filter_for_zilany2009(sound, fs)



    psp = []
    for freq in cfs:

        # Run IHC model
        vihc = _pycat.run_ihc(
            signal=sound,
            cf=float(freq),
            fs=float(fs),
            cohc=float(cohc),
            cihc=float(cihc)
        )


        synout = _pycat.run_synapse(
            fs=fs,
            vihc=vihc,
            cf=freq,
            anf_type=anf_type,
            powerlaw_implnt=powerlaw_implnt,
            with_ffGn=with_ffGn
        )

        psp.append(synout)

    psp = np.array(psp).T


    return psp
