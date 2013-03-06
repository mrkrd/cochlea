#!/usr/bin/env python

"""Holmberg, M. (2007). Speech Encoding in the Human Auditory
Periphery: Modeling and Quantitative Assessment by Means of Automatic
Speech Recognition. PhD thesis, Technical University Darmstadt.

"""


from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import os
import pandas as pd
import itertools

import traveling_waves as tw


def par_dir(par_file):
    """
    Add directory path to par file name that is refered to the
    modules location.
    """
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'pars', par_file)


def run_holmberg2007(
        sound,
        fs,
        anf_num,
        seed,
        cf=None
):

    assert np.max(np.abs(sound)) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1
    assert cf is None
    assert fs == 48e3

    np.random.seed(seed)



    duration = len(sound) / fs



    ### Outer ear filter
    sound_oe = tw.run_outer_ear_filter(sound, fs)

    ### Middle ear filter
    sound_me = tw.run_middle_ear_filter(sound_oe, fs)

    ### Scaling
    sound_scaled = sound_me * tw.scaling_factor

    ### Basilar membrane
    xbms = tw.run_bm_wave(sound_scaled, fs)


    ihcrp = {}
    for cf,xbm in xbms.iteritems():

        ### Amplification
        lcr4 = tw.run_lcr4(xbm, fs, cf)

        ### IHCRP
        ihcrp[cf] = tw.run_ihcrp(lcr4, fs, cf)



    anf_types = np.repeat(['hsr', 'msr', 'lsr'], anf_num)

    psps = {}
    trains = []
    for cf,anf_type in itertools.product(ihcrp.keys(),anf_types):

        if (cf,anf_type) not in psps:
            psp = tw.run_ihc_meddis2000(
                ihcrp=ihcrp[cf],
                fs=fs,
                gamma_Ca=130,
                beta_Ca=400,
                tau_m=1e-4,
                G_Ca_max=4.5e-9,
                E_Ca=0.066,
                tau_Ca=1e-4,
                perm_Ca0=0,
                perm_z=2e32,
                pCa=3,
                replenish_rate_y=10,
                loss_rate_l=2580,
                recovery_rate_r=6580,
                reprocess_rate_x=66.3,
                max_free_pool=8,
            )
            psps[cf,anf_type] = psp

        spikes = tw.run_an_sg_carney_holmberg2007(
            psp=psp,
            fs=fs,
            c0=0.5,
            c1=0.5,
            s0=1e-3,
            s1=12.5e-3,
            refractory_period=0.75e-3
        )
        trains.append({
            'spikes': spikes,
            'duration': duration,
            'cf': cf,
            'type': anf_type
        })


    spike_trains = pd.DataFrame(trains)
    return spike_trains
