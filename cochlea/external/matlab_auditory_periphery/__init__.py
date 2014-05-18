#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np

import pymatlab

import mrlib.thorns as th


def run_matlab_auditory_periphery(
        sound,
        fs,
        anf_num,
        cf,
        seed,
):
    """Wrapper function for an auditory model implemented in Matlab
Auditory Periphery (Meddis et al.).

    Requires Matlab and pymatlab.

    """

    ### Validate `anf_num`
    assert len(anf_num) == 3
    if not (anf_num[0] == anf_num[1] == anf_num[2]):
        raise RuntimeError("All elements of anf_num must be equal: {}".format(anf_num))


    ### Validate `cf`
    if isinstance(cf, tuple):
        assert len(cf) == 3
    elif len(cf) == 3:
        raise RuntimeError("Three frequency channels are forbidden, because they mask the tuple (min_cf, max_cf, cf_num).")
    elif np.isscalar(cf):
        cf = [float(cf)]

    matlab = pymatlab.session_factory('-nojvm')

    matlab.run("rng({})".format(seed))

    matlab.run("global dtSpikes ANoutput savedBFlist")

    matlab.putvalue("inputSignal", sound)
    matlab.putvalue("sampleRate", np.array([fs], dtype=float))
    matlab.putvalue("BFlist", np.array(cf, dtype=float))


    cmd = "MAP1_14(inputSignal, sampleRate, BFlist, 'Normal', 'spikes', {{'AN_IHCsynapseParams.numFibers={anf_num};','AN_IHCsynapseParams.spikesTargetSampleRate={fs};'}})".format(
        anf_num=anf_num[0],
        fs=fs,
    )

    matlab.run(cmd)

    anf_raw = matlab.getvalue("ANoutput")
    cf_raw = matlab.getvalue("savedBFlist")
    dt_raw = matlab.getvalue("dtSpikes")


    anf_trains = th.make_trains(
        anf_raw.T,
        fs=1/dt_raw
    )

    cf_col = np.tile(np.repeat(cf_raw,anf_num[0]), 3)
    type_col = np.repeat(['lsr', 'msr', 'hsr'], len(cf_raw)*anf_num[0])

    anf_trains['type'] = type_col
    anf_trains['cf'] = cf_col


    return anf_trains






def main():
    pass

if __name__ == "__main__":
    main()
