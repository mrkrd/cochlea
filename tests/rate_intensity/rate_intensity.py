# Author: Marek Rudnicki
# Time-stamp: <2009-12-17 18:24:45 marek>
#
# Description: Rate-intensity function

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import cochlea
import thorns as th

def rate_intensity(ear, fs, cf_list, dbspl_list):

    tmax = 0.1
    t = np.arange(0, tmax, 1/fs)

    hsr_rate = np.zeros( (len(cf_list), len(dbspl_list)) )
    msr_rate = np.zeros( (len(cf_list), len(dbspl_list)) )
    lsr_rate = np.zeros( (len(cf_list), len(dbspl_list)) )


    for cf_idx,cf in enumerate(cf_list):
        s1 = np.sin(2 * np.pi * t * cf)
        ear.set_freq(cf)

        for dbspl_idx,dbspl in enumerate(dbspl_list):
            s = cochlea.set_dB_SPL(dbspl, s1)

            print cf, dbspl
            hsr, msr, lsr = ear.run(fs, s)

            if hsr != None:
                all_spikes = np.concatenate(tuple(hsr['spikes']))
                hsr_rate[cf_idx,dbspl_idx] = len(all_spikes) / tmax / ear.hsr_num
            if msr != None:
                all_spikes = np.concatenate(tuple(msr['spikes']))
                msr_rate[cf_idx,dbspl_idx] = len(all_spikes) / tmax / ear.msr_num
            if lsr != None:
                all_spikes = np.concatenate(tuple(lsr['spikes']))
                lsr_rate[cf_idx,dbspl_idx] = len(all_spikes) / tmax / ear.lsr_num


    return hsr_rate, msr_rate, lsr_rate




def rate_intensity_sumner2002():

    dbspl_list=range(-20,80,5)
    cf_list = [10000]

    ear = cochlea.Sumner2002(hsr=10, msr=10, lsr=10)

    hsr_rate, msr_rate, lsr_rate = rate_intensity(ear, 100000,
                                                  cf_list=cf_list,
                                                  dbspl_list=dbspl_list)
    np.save('hsr_rate.npy', hsr_rate)
    plt.plot(dbspl_list, hsr_rate[0])
    plt.plot(dbspl_list, msr_rate[0])
    plt.plot(dbspl_list, lsr_rate[0])
    plt.show()

def rate_intensity_carney2009():

    dbspl_list=range(-20,80,5)
    cf_list = [1000]

    ear = cochlea.Carney2009(hsr=50, msr=50, lsr=50,
                             powerlaw_implnt='approx')

    hsr_rate, msr_rate, lsr_rate = rate_intensity(ear, 100000,
                                                  cf_list=cf_list,
                                                  dbspl_list=dbspl_list)
    np.save('hsr_rate.npy', hsr_rate)

    fig = plt.gcf()
    ax = fig.add_subplot(111)

    ax.plot(dbspl_list, hsr_rate[0])
    ax.plot(dbspl_list, msr_rate[0])
    ax.plot(dbspl_list, lsr_rate[0])

    ax.set_xlabel("Intensity (dB SPL)")
    ax.set_ylabel("Rate (spikes / sec)")

    fig.savefig('carney2009-rate-intensity.eps')

    plt.show()


if __name__ == "__main__":
    rate_intensity_carney2009()
