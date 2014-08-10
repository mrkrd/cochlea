# -*- coding: utf-8 -*-

# Copyright 2014 Marek Rudnicki
#
# This file is part of cochlea.
#
# cochlea is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cochlea is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cochlea.  If not, see <http://www.gnu.org/licenses/>.


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
        params_name='Normal',
        matlab_session=None
):
    """Run Matlab Auditory Periphery [MAP]_ model by Ray Meddis.  This
    function does not implement the model, but wraps the model
    implementation using `matlab_wrapper`.  The model takes sound
    signal as input and outputs auditory nerve spike trains.

    In order to run it, make sure that all necessary [MAP]_ model
    files are in `MATLABPATH`.  You should be able to run `MAP1_14`
    function in MATLAB first!

    Requires MATLAB, *matlab_wrapper*, *thorns* and [MAP]_ in
    `MATLABPATH`.


    Parameters
    ----------
    sound : array_like
        Input sound.
    fs : float
        Sampling frequency of the sound.
    anf_num : tuple
        Number of auditory nerve fibers per channel (HSR#, MSR#, LSR#).
    cf : float or array_like or tuple
        The center frequency(s) of the simulated auditory nerve fibers.
        If float, then defines a single frequency channel.  If
        array_like (e.g. list or ndarray), then the frequencies are
        used.  If tuple, then must have exactly 3 elements (min_cf,
        max_cf, num_cf) and the frequencies are calculated using the
        Greenwood function.
    seed : int
        Random seed.
    params_name : str, optional Tail of the parameter filename
        (parameterStore/MAPparams<params_name>.m).  Refer to MAP
        documentation for the details.
    matlab_session : MatlabSession or None, optional
        MatlabSession object from `matlab_wrapper` module.  If `None`,
        then new session is generated.


    Returns
    -------
    spike_trains
        Auditory nerve spike trains.


    References
    ----------

    .. [MAP] http://www.essexpsychology.macmate.me/HearingLab/modelling.html

    """

    ### Validate `anf_num`
    assert len(anf_num) == 3


    ### Validate `cf`
    if isinstance(cf, tuple):
        assert len(cf) == 3
    elif np.isscalar(cf):
        pass
    elif len(cf) == 3:
        raise RuntimeError("Three frequency channels are forbidden, because they mask the tuple (min_cf, max_cf, cf_num).")


    ### Generate matlab_wrapper session as needed
    if matlab_session is None:
        matlab = matlab_wrapper.MatlabSession(options='-nosplash -singleCompThread')
    else:
        matlab = matlab_session

    matlab.eval("clear all")
    matlab.eval("clearvars -global")

    ### Set Matlab environment
    matlab.workspace.rng(seed)

    matlab.eval("global dtSpikes ANoutput savedBFlist")

    matlab.workspace.MAP1_14(
        sound,
        float(fs),
        np.array(cf, dtype=float),
        params_name,
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
