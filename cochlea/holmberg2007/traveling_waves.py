from __future__ import division

import numpy as np
import scipy.signal as dsp

from cochlea.holmberg2007.bm_pars import (
    real_freq_map,
    S_ST,
    S_ED,
    C_eardrum,
    outer_ear_a_48kHz,
    outer_ear_b_48kHz,
    outer_ear_a_100kHz,
    outer_ear_b_100kHz,
)
from _traveling_waves import *


# Input signal should be multiplied by this factor
scaling_factor = S_ST * S_ED



def run_middle_ear_filter_orig(signal, fs):
    """ Middle ear filter designed using digital wave techinique. """
    assert fs == 48000

    R2_ME = 1. / (2. * fs * C_eardrum)
    R1 = 1. / (2. * np.pi * C_eardrum * 1e3)

    g1_ME = (R2_ME - R1) / (R2_ME + R1)
    Z_ME=0

    out = np.zeros(len(signal))

    for i,samp in enumerate(signal):
        b2 = samp + g1_ME * (samp - Z_ME)
        Z_ME_b = Z_ME
        Z_ME = b2
        out[i] = ((Z_ME_b - b2) / R2_ME)

    return out


def _calc_middle_ear_coefs(fs):
    # Digital wave filter coefficients
    R2_ME = 1. / (2. * fs * C_eardrum)
    R1 = 1. / (2. * np.pi * C_eardrum * 1e3)
    g1_ME = (R2_ME - R1) / (R2_ME + R1)
    Z_ME=0

    # Standard filter coefficients
    b = [(-1-g1_ME), (1+g1_ME)]
    a = [R2_ME, R2_ME*g1_ME]

    return b, a


def run_middle_ear_filter(signal, fs):
    """ Middle ear filter model. """
    b,a = _calc_middle_ear_coefs(fs)

    return dsp.lfilter(b, a, signal)


def _calc_outer_ear_coefs(fs):
    if fs == 48000:
        a = outer_ear_a_48kHz
        b = outer_ear_b_48kHz
    elif fs == 100000:
        a = outer_ear_a_100kHz
        b = outer_ear_b_100kHz
    else:
        assert False, "Invalid sampling frequency: fs"

    return b, a


def run_outer_ear_filter(signal, fs):
    b, a = _calc_outer_ear_coefs(fs)

    return dsp.lfilter(b, a, signal)



def find_closest_freq_idx_in_map(freq):
    m = np.abs(real_freq_map - freq)

    return int(np.argmin(m))


def find_closest_freq_in_map(freq):
    idx = find_closest_freq_idx_in_map(freq)
    return real_freq_map[idx]



def run_ihc_meddis2000(
        ihcrp,
        fs,
        gammaCa,
        betaCa,
        tauCaChan,
        GCaMax,
        CaVrev,
        tauConcCa,
        perm_Ca0,
        perm_z,
        pCa,
        replenishRate_y,
        lossRate_l,
        recoveryRate_r,
        reprocessRate_x,
        maxFreePool_M,
        opmode='probability'
):

    dt = 1/fs

    psp = np.zeros_like(ihcrp)

    ### Initial values
    uIHC_resting = ihcrp[0]


    ssactCa = 1 / (1 + np.exp(-gammaCa*uIHC_resting)/betaCa)
    ICa = GCaMax * ssactCa**3 * (uIHC_resting - CaVrev)

    if -ICa > perm_Ca0:
        spontPerm_k0 = perm_z * ((-ICa)**pCa - perm_Ca0**pCa)
    else:
        spontPerm_k0 = 0.0

    # cleftReplenishMode == IHC_MEDDIS2000_CLEFTREPLENISHMODE_ORIGINAL
    spontCleft_c0 = (
        maxFreePool_M * replenishRate_y * spontPerm_k0
        /
        (replenishRate_y * (lossRate_l+recoveryRate_r) + spontPerm_k0*lossRate_l)
    )

    if (spontCleft_c0 > 0) and opmode == 'probability':
        spontFreePool_q0 = spontCleft_c0 * (lossRate_l + recoveryRate_r) / spontPerm_k0
    elif (spontCleft_c0 > 0) and opmode == 'spikes':
        spontFreePool_q0 = np.floor( (spontCleft_c0 * (lossRate_l+recoveryRate_r) / spontPerm_k0) + 0.5 )
    else:
        spontFreePool_q0 = maxFreePool_M

    spontReprocess_w0 = spontCleft_c0 * recoveryRate_r / reprocessRate_x

    actCa = ssactCa
    concCa = -ICa
    reservoirQ = spontFreePool_q0
    cleftC = spontCleft_c0
    reprocessedW = spontReprocess_w0



    for i,Vin in enumerate(ihcrp):


        ### Ca current
        ssactCa = 1 / (1 + np.exp(-gammaCa*Vin)/betaCa)
        actCa += (ssactCa - actCa) * dt / tauCaChan
        ICa = GCaMax * actCa**3 * (Vin - CaVrev)


        ### Calcium Ion accumulation and diffusion
        # caCondMode == IHC_MEDDIS2000_CACONDMODE_ORIGINAL
        concCa += (-ICa - concCa) * dt / tauConcCa


        ### power law release function
        if concCa > perm_Ca0:
            kdt = (perm_z*dt * (concCa**pCa - perm_Ca0**pCa))
        else:
            kdt = 0.0



        ### Synapse
        if opmode == 'probability':
            # cleftReplenishMode == IHC_MEDDIS2000_CLEFTREPLENISHMODE_ORIGINAL
            if reservoirQ < maxFreePool_M:
                replenish = replenishRate_y * dt * (maxFreePool_M - reservoirQ)
            else:
                replenish = 0.0


            ejected = kdt * reservoirQ

            reUptakeAndLost = (lossRate_l + recoveryRate_r) * dt * cleftC
            reUptake = recoveryRate_r * dt * cleftC

            reprocessed = reprocessRate_x * dt * reprocessedW
            reservoirQ += replenish - ejected + reprocessed
            cleftC += ejected - reUptakeAndLost

            psp[i] = ejected

            reprocessedW += reUptake - reprocessed
        elif opmode == 'spikes':
            if reservoirQ < maxFreePool_M:
                replenish = (np.random.geometric(
                    replenishRate_y * dt,
                    int(maxFreePool_M - reservoirQ)
                ) == 1).sum()
            else:
                replenish = 0.0

            ejected = (np.random.geometric(
                kdt,
                int(reservoirQ)
            ) == 1).sum()

            reUptakeAndLost = (lossRate_l+recoveryRate_r)*dt * cleftC
            reUptake = recoveryRate_r * dt * cleftC
            if reprocessedW < 1:
                reprocessed = 0.0
            else:
                reprocessed = (np.random.geometric(
                    reprocessRate_x * dt,
                    int(reprocessedW)
                ) == 1).sum()

            reservoirQ += replenish - ejected + reprocessed
            cleftC += ejected - reUptakeAndLost


            if ejected > 0:
                psp[i] = ejected;
            else:
                psp[i] = 0.0

            reprocessedW += reUptake - reprocessed


        else:
            raise RuntimeError

    return psp
