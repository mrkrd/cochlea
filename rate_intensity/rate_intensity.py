# Author: Marek Rudnicki
# Time-stamp: <2009-09-15 21:39:08 marek>
#
# Description: Rate-intensity function


import numpy as np
import matplotlib.pyplot as plt

import cochlea

def rate_intensity(ear, fs, cf_list, dbspl_list):

    fs = float(fs)

    tmax = 0.1
    t = np.arange(0, tmax, 1./fs)

    hsr_rate = np.zeros( (len(cf_list), len(dbspl_list)) )
    msr_rate = np.zeros( (len(cf_list), len(dbspl_list)) )
    lsr_rate = np.zeros( (len(cf_list), len(dbspl_list)) )


    for cf_idx,cf in enumerate(cf_list):
        s1 = np.sin(2 * np.pi * t * cf)
        ear.set_freq(cf)

        for dbspl_idx,dbspl in enumerate(dbspl_list):
            s = cochlea.set_dB_SPL(dbspl, s1)

            hsr, msr, lsr = ear.run(fs, s)

            if hsr != None:
                all_spikes = np.concatenate(tuple(hsr['spikes']))
                hsr_rate[cf_idx,dbspl_idx] = len(all_spikes) / tmax
            if msr != None:
                all_spikes = np.concatenate(tuple(msr['spikes']))
                msr_rate[cf_idx,dbspl_idx] = len(all_spikes) / tmax
            if lsr != None:
                all_spikes = np.concatenate(tuple(lsr['spikes']))
                lsr_rate[cf_idx,dbspl_idx] = len(all_spikes) / tmax


    return hsr_rate, msr_rate, lsr_rate




def rate_intensity_sumner2002():

    dbspl_list=range(0,80,2)

    ear = cochlea.Sumner2002(hsr=20, msr=0, lsr=0)

    hsr_rate, msr_rate, lsr_rate = rate_intensity(ear, 100000,
                                                  cf_list=[10000, 30000],
                                                  dbspl_list=dbspl_list)
    np.save('hsr_rate.npy', hsr_rate)
    plt.plot(dbspl_list, hsr_rate.T)
    plt.show()


if __name__ == "__main__":
    rate_intensity_sumner2002()
