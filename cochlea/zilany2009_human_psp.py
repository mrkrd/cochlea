# Description: Model of auditory periphery of: Zilany, M.S.A., Bruce,
# I.C., Nelson, P.C., and Carney, L.H. (manuscript in preparation) 2009

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import scipy.signal as dsp
import os

from cochlea.pycat import _pycat
from cochlea.zilany2009_human import (
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
        parallel=False):



    assert np.max(sound) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1


    cfs = _calc_cfs(cf)


    # Run Outer/Middle Ear filter
    sound = _run_human_me_filter_for_zilany2009(sound, fs)


    channel_args = [
        {
            'sound': sound,
            'cf': freq,
            'fs': fs,
            'cohc': cohc,
            'cihc': cihc,
            'anf_type': anf_type,
            'powerlaw_implnt': powerlaw_implnt,
        }
        for freq in cfs
    ]


    if parallel:
        import multiprocessing

        pool = multiprocessing.Pool()
        psp = pool.map(_run_channel, channel_args)

    else:
        psp = map(_run_channel, channel_args)


    psp = np.array(psp).T


    return psp, cfs



def _run_channel(args):

    vihc = _pycat.run_ihc(
        signal=args['sound'],
        cf=float(args['cf']),
        fs=float(args['fs']),
        cohc=float(args['cohc']),
        cihc=float(args['cihc'])
    )


    synout = _pycat.run_synapse(
        fs=args['fs'],
        vihc=vihc,
        cf=args['cf'],
        anf_type=args['anf_type'],
        powerlaw_implnt=args['powerlaw_implnt'],
        with_ffGn=False
    )

    return synout
