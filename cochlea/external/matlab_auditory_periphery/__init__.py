#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np
import pandas as pd

import matlab_wrapper

import thorns as th


def run_matlab_auditory_periphery(
        sound,
        fs,
        anf_num,
        cf,
        seed,
        matlab_session=None
):
    """Wrapper function for an auditory model implemented in Matlab
    Auditory Periphery (Meddis et al.).

    Requires Matlab and matlab_wrapper.

    """

    ### Validate `anf_num`
    assert len(anf_num) == 3


    ### Validate `cf`
    if isinstance(cf, tuple):
        assert len(cf) == 3
    elif len(cf) == 3:
        raise RuntimeError("Three frequency channels are forbidden, because they mask the tuple (min_cf, max_cf, cf_num).")


    ### Generate matlab_wrapper session as needed
    if matlab_session is None:
        matlab = matlab_wrapper.MatlabSession(options='-nosplash -singleCompThread', buffer_size=1000)
    else:
        matlab = matlab_session


    ### Set Matlab environment
    matlab.workspace.rng(seed)

    matlab.eval("global dtSpikes ANoutput savedBFlist")

    matlab.workspace.MAP1_14(
        sound,
        float(fs),
        np.array(cf, dtype=float),
        'Normal',
        'spikes',
        ['AN_IHCsynapseParams.numFibers={};'.format(max(anf_num)),
         'AN_IHCsynapseParams.spikesTargetSampleRate={};'.format(fs)],
        nout=0
    )


    ### Collect results from Matlab
    anf_raw = matlab.workspace.ANoutput

    cf_raw = matlab.workspace.savedBFlist
    cf_raw = np.atleast_1d(cf_raw)

    dt_raw = matlab.workspace.dtSpikes



    ### Make trains
    all_trains = th.make_trains(
        anf_raw.T,
        fs=1/dt_raw
    )

    cf_col = np.tile(np.repeat(cf_raw,anf_num[0]), 3)
    type_col = np.repeat(['lsr', 'msr', 'hsr'], len(cf_raw)*anf_num[0])

    all_trains['type'] = type_col
    all_trains['cf'] = cf_col



    n = {
        'hsr': anf_num[0],
        'msr': anf_num[1],
        'lsr': anf_num[2],
    }

    ### Discard trains.  We want only anf_num == (HSR#, MSR#, LSR#)
    anf_trains = []
    for name,group in all_trains.groupby(['type','cf']):
        typ,cf = name

        sel = group.iloc[0:n[typ]]

        anf_trains.append(sel)


    anf_trains = pd.concat(anf_trains)

    return anf_trains
