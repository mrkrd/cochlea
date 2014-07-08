"""Copyright 2009-2014 Marek Rudnicki

This file is part of cochlea.

cochlea is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

cochlea is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with cochlea.  If not, see <http://www.gnu.org/licenses/>.

"""

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"


import numpy as np
import pandas as pd
import itertools

from . import traveling_waves as tw
from . traveling_waves import real_freq_map, get_nearest_cf


def run_holmberg2007(
        sound,
        fs,
        anf_num,
        seed,
        cf=None,
):
    """Run the inner ear model by [Holmberg2007]_.  It simulates the
    traveling wave on the basilar membrane, inner hair cell, synapses
    and generates auditory nerve spikes.  The model takes sound signal
    as input and outputs auditory nerve spike trains.


    Parameters
    ----------
    sound : array_like
        Input sound signal.
    fs : float
        Sampling frequency of the sound.
    anf_num : tuple
        Number of auditory nerve fibers per channel (HSR#, MSR#, LSR#).
    seed : int
        Random seed.
    cf : float or array_like or None, optional
        Characteristic frequencies.  If `None`, then calculate all 100
        predefined CFs.  CFs must be a subset of `real_freq_map`.


    Returns
    -------
    spike_trains
        Auditory nerve spike trains.


    References
    ----------
    If you are using results of this or modified version of the model
    in your research, please cite [Holmberg2007]_.

    .. [Holmberg2007] Holmberg, M. (2007). Speech Encoding in the
       Human Auditory Periphery: Modeling and Quantitative Assessment
       by Means of Automatic Speech Recognition. PhD thesis, Technical
       University Darmstadt.

    """

    assert np.max(np.abs(sound)) < 1000, "Signal should be given in Pa"
    assert sound.ndim == 1
    assert fs == 48e3

    np.random.seed(seed)

    if cf is None:
        cfs = tw.real_freq_map
    elif np.isscalar(cf):
        cfs = [cf]
    else:
        cfs = cf


    assert set(cfs) <= set(tw.real_freq_map), set(cfs) - set(tw.real_freq_map)


    duration = len(sound) / fs



    ### Outer ear filter
    sound_oe = tw.run_outer_ear_filter(sound, fs)

    ### Middle ear filter
    sound_me = tw.run_middle_ear_filter(sound_oe, fs)

    ### Scaling
    sound_scaled = sound_me * tw.scaling_factor

    ### Basilar membrane
    xbm = tw.run_bm_wave(sound_scaled, fs)


    ihcrp = {}
    for cf in cfs:
        ### Amplification
        lcr4 = tw.run_lcr4(xbm[cf], fs, cf)

        ### Delay correction (1/cf)
        sec = np.where(bm_pars.real_freq_map == cf)[0][0]
        lcr4_rolled = np.roll(
            lcr4,
            -int(np.round(tw.delay_time[99-sec]*fs))
        )

        ### IHCRP
        ihcrp[cf] = tw.run_ihcrp(lcr4_rolled, fs, cf)



    anf_types = np.repeat(['hsr', 'msr', 'lsr'], anf_num)


    ihc_meddis2000_pars = {
        'hsr': {
            'gamma_ca': 130,
            'beta_ca': 400,
            'tau_m': 1e-4,
            'gmax_ca': 4.5e-9,
            'e_ca': 0.066,
            'tau_ca': 1e-4,
            'perm_ca0': 0,
            'perm_z': 2e32,
            'power_ca': 3,
            'replenish_rate_y': 10,
            'loss_rate_l': 2580,
            'recovery_rate_r': 6580,
            'reprocess_rate_x': 66.3,
            'max_free_pool': 8,
        },
        'msr': {
            'gamma_ca': 130,
            'beta_ca': 400,
            'tau_m': 1e-4,
            'gmax_ca': 4.25e-9,
            'e_ca': 0.066,
            'tau_ca': 1e-4,
            'perm_ca0': 2.5e-11,
            'perm_z': 2e32,
            'power_ca': 3,
            'replenish_rate_y': 10,
            'loss_rate_l': 2580,
            'recovery_rate_r': 6580,
            'reprocess_rate_x': 66.3,
            'max_free_pool': 9,
        },
        'lsr': {
            'gamma_ca': 130,
            'beta_ca': 400,
            'tau_m': 1e-4,
            'gmax_ca': 2.75e-9,
            'e_ca': 0.066,
            'tau_ca': 1e-4,
            'perm_ca0': 4.2e-11,
            'perm_z': 2e32,
            'power_ca': 3,
            'replenish_rate_y': 10,
            'loss_rate_l': 2580,
            'recovery_rate_r': 6580,
            'reprocess_rate_x': 66.3,
            'max_free_pool': 6,
        }
    }

    psps = {}
    trains = []
    for cf,anf_type in itertools.product(ihcrp.keys(),anf_types):

        if (cf,anf_type) not in psps:
            ### IHC and Synapse
            psp = tw.run_ihc_meddis2000(
                ihcrp=ihcrp[cf],
                fs=fs,
                **ihc_meddis2000_pars[anf_type]
            )
            psps[cf,anf_type] = psp

        ### Spike generator (pars from Sumner et al. 2002)
        spikes = tw.run_an_sg_carney_holmberg2007(
            psp=psp,
            fs=fs,
            c0=0.5,
            c1=0.,
            s0=0.8e-3,
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
