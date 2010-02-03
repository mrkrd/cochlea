# Author: Marek Rudnicki
# Time-stamp: <2010-02-03 14:32:47 marek>
#
# Description: Rate-intensity function

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import cochlea
import thorns as th
import thorns.waves as wv

def rate_intensity(ear, fs, cf, dbspl_list):

    tmax = 100
    s = wv.generate_ramped_tone(fs, cf,
                                tone_duration=tmax,
                                pad_duration=0)
    ear.set_freq(cf)

    hsr_rates = []
    msr_rates = []
    lsr_rates = []

    for dbspl in dbspl_list:
        s = wv.set_dbspl(dbspl, s)

        print "dB_SPL:", dbspl
        hsr, msr, lsr = ear.run(fs, s)

        if hsr != None:
            hsr_rates.append( th.firing_rate( hsr['spikes'],
                                              tmax, ear._hsr_num) )
        if msr != None:
            msr_rates.append( th.firing_rate( msr['spikes'],
                                              tmax, ear._msr_num) )
        if lsr != None:
            lsr_rates.append( th.firing_rate( lsr['spikes'],
                                              tmax, ear._lsr_num) )

    return hsr_rates, msr_rates, lsr_rates




def rate_intensity_sumner2002():
    dbspl_list = np.arange(-20, 120, 5)
    cf = 5000

    ear = cochlea.Sumner2002((50, 50, 50), freq=cf)

    hsr_rates, msr_rates, lsr_rates = \
        rate_intensity(ear, fs=100000, cf=cf, dbspl_list=dbspl_list)

    fig = plt.gcf()
    ax = fig.add_subplot(111)

    ax.plot(dbspl_list, hsr_rates)
    ax.plot(dbspl_list, msr_rates)
    ax.plot(dbspl_list, lsr_rates)

    ax.set_xlabel("Intensity (dB SPL)")
    ax.set_ylabel("Rate (spikes / sec)")

    fig.savefig("sumner2002_rate-intensity.eps")

    plt.show()


def rate_intensity_holmberg2008():
    import traveling_waves as tw

    dbspl_list = np.arange(-20, 120, 5)
    cf = tw.real_freq_map[50]

    ear = cochlea.Holmberg2008((50, 50, 50), freq=cf)

    hsr_rates, msr_rates, lsr_rates = \
        rate_intensity(ear, fs=48000, cf=cf, dbspl_list=dbspl_list)

    fig = plt.gcf()
    ax = fig.add_subplot(111)

    ax.plot(dbspl_list, hsr_rates)
    ax.plot(dbspl_list, msr_rates)
    ax.plot(dbspl_list, lsr_rates)

    ax.set_xlabel("Intensity (dB SPL)")
    ax.set_ylabel("Rate (spikes / sec)")

    fig.savefig("holmberg2008_rate-intensity.eps")

    plt.show()



def rate_intensity_zilany2009():

    dbspl_list = np.arange(-20, 120, 5)
    cf = 2000

    ear = cochlea.Zilany2009((50, 50, 50), powerlaw_implnt='approx')

    hsr_rates, msr_rates, lsr_rates = \
        rate_intensity(ear, fs=100000, cf=cf, dbspl_list=dbspl_list)

    fig = plt.gcf()
    ax = fig.add_subplot(111)

    ax.plot(dbspl_list, hsr_rates)
    ax.plot(dbspl_list, msr_rates)
    ax.plot(dbspl_list, lsr_rates)

    ax.set_xlabel("Intensity (dB SPL)")
    ax.set_ylabel("Rate (spikes / sec)")

    fig.savefig("zilany2009_rate-intensity.eps")

    plt.show()


if __name__ == "__main__":
    rate_intensity_zilany2009()
    # rate_intensity_sumner2002()
    # rate_intensity_holmberg2008()
