#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt

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


def main():

    fs = 48e3

    ihcrp = np.loadtxt(
        "data/ihcrp_fs_48e3.csv",
    )

    synout = run_ihc_meddis2000(
        ihcrp=ihcrp,
        fs=fs,
        gammaCa=130,
        betaCa=400,
        tauCaChan=1e-4,
        GCaMax=4.5e-9,
        CaVrev=0.066,
        tauConcCa=1e-4,
        perm_Ca0=0,
        perm_z=2e32,
        pCa=3,                  # default
        replenishRate_y=10,
        lossRate_l=2580,
        recoveryRate_r=6580,
        reprocessRate_x=66.3,
        maxFreePool_M=8,
        # opmode='spikes'
    )

    target = np.loadtxt(
        "data/ihc_meddis2000_dsam.csv",
    )

    plt.plot(synout)
    plt.plot(target)
    plt.show()


if __name__ == "__main__":
    main()
