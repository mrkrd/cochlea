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

import traveling_waves as tw
import dsam
import marlib.thorns as th

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

    assert np.max(sound) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1
    assert cf is None
    assert fs == 48e3

    np.random.seed(seed)




    freq_map = tw.bm_pars.real_freq_map
    duration = len(sound) / fs





    ### Outer ear filter
    sound_oe = tw.run_outer_ear_filter(sound, fs)

    ### Middle ear filter
    sound_me = tw.run_middle_ear_filter(sound_oe, fs)

    ### Scaling
    sound_scaled = sound_me * tw.scaling_factor

    ### Basilar membrane
    xbm = tw.run_bm_wave(sound_scaled, fs)

    ### Amplification
    LCR4 = tw.run_LCR4(xbm, fs)

    ### IHCRP
    ihcrp = tw.run_ihcrp(LCR4, fs)


    anf_types = np.repeat(['hsr', 'msr', 'lsr'], anf_num)
    ihcs = {'hsr':None, 'msr':None, 'lsr':None}
    sgs = {'hsr':None, 'msr':None, 'lsr':None}
    trains = []

    for anf_type in anf_types:

        if ihcs[anf_type] is None:
            ihc_module = dsam.EarModule("IHC_Meddis2000")
            ihc_module.read_pars(par_dir("ihc_%s_Sumner2002.par" %anf_type))
            sg_module = dsam.EarModule("An_SG_Carney")
            sg_module.read_pars(par_dir("anf_carney.par"))
            sg_module.set_par('pulse_duration', 1.1/fs)
            dsam.connect(ihc_module, sg_module)

            ihc_module.run(fs, ihcrp)

            ihcs[anf_type] = ihc_module
            sgs[anf_type] = sg_module


        sgs[anf_type].run()
        sg_signal = sgs[anf_type].get_signal()


        ### Convert array to `nude' spike trains (no meta)
        nudes = []
        for sig in sg_signal.T:
            sig = sig.astype(int)
            t = np.arange(len(sig))
            spikes = np.repeat(t, sig) / fs

            nudes.append(spikes)



        for cf,spikes in zip(freq_map, nudes):
            trains.append({
                'spikes': spikes,
                'duration': duration,
                'cf': cf,
                'type': anf_type
            })



    spike_trains = pd.DataFrame(trains)
    return spike_trains
