# Author: Marek Rudnicki
# Time-stamp: <2009-09-09 14:10:24 marek>
#
# Description: Testing various SG modules

import numpy as np
import matplotlib.pyplot as plt

import dsam
import cochlea
import thorns as th

def main():
    fs = 48000.
    cf = 4000.

    nfib = 1000

    t = np.arange(0, 0.1, 1./fs)
    s = np.sin(2 * np.pi * t * cf)
    s = cochlea.set_dB_SPL(60, s)

    ear = cochlea.Sumner2002(hsr=1, msr=0, lsr=0,
                             freq = cf,
                             animal='human')

    ear.run(fs, s)

    spike_prob = ear.ihc_hsr.get_signal()



    sg_simple = dsam.EarModule("AN_SG_Simple")
    sg_simple.set_par("NUM_FIBRES", nfib)
    sg_simple.set_par("PULSE_DURATION", 1.1/fs)
    sg_simple.set_par("REFRAC_PERIOD", 0.00075)
    sg_simple.set_par("MAGNITUDE", 1)
    sg_simple.print_pars()
    sg_simple.run(fs, spike_prob)


    sg_binomial = dsam.EarModule("AN_SG_Binomial")
    sg_binomial.set_par("NUM_FIBRES", nfib)
    sg_binomial.set_par("PULSE_DURATION", 1.1/fs)
    sg_binomial.set_par("REFRAC_PERIOD", 0.00075)
    sg_binomial.set_par("MAGNITUDE", 1)
    sg_binomial.print_pars()
    sg_binomial.run(fs, spike_prob)

    # sg_meddis = dsam.EarModule("AN_SG_Meddis02")
    # sg_meddis.set_par("NUM_FIBRES", nfib)
    # sg_meddis.set_par("PULSE_DURATION", 1.1/fs)
    # sg_meddis.set_par("REFRAC_PERIOD", 0.00075)
    # sg_meddis.set_par('RECOVERY_TAU', 0.0008)
    # sg_meddis.set_par("MAGNITUDE", 1)
    # sg_meddis.print_pars()
    # sg_meddis.run(fs, spike_prob)

    sg_carney = dsam.EarModule("AN_SG_Carney")
    sg_carney.set_par("NUM_FIBRES", nfib)
    #sg_carney.set_par("PULSE_DURATION", 1./fs)
    sg_carney.set_par("REFRAC_PERIOD", 0.00075)
    sg_carney.set_par("C0", 0.55)
    sg_carney.set_par("C1", 0.0)
    sg_carney.set_par("S0", 0.0008)
    sg_carney.set_par("S1", 0.0008)
    sg_carney.set_par("MAGNITUDE", 1)
    sg_carney.print_pars()
    sg_carney.run(fs, spike_prob)

    simple_spikes = th.signal_to_spikes(fs, sg_simple.get_signal())
    binomial_spikes = th.signal_to_spikes(fs, sg_binomial.get_signal())
    # meddis_spikes = th.signal_to_spikes(fs, sg_meddis.get_signal())
    carney_spikes = th.signal_to_spikes(fs, sg_carney.get_signal())

    print th.synchronization_index(cf, simple_spikes)
    print th.synchronization_index(cf, binomial_spikes)
    print th.synchronization_index(cf, carney_spikes)


    ax = plt.gca()

    th.plot_psth(simple_spikes, bin_size=0.1, axis=ax, histtype='step')
    th.plot_psth(binomial_spikes, bin_size=0.1, axis=ax, histtype='step')
    # th.plot_psth(meddis_spikes, bin_size=0.1, axis=ax, histtype='step')
    th.plot_psth(carney_spikes, bin_size=0.1, axis=ax, histtype='step')

    plt.show()


if __name__ == "__main__":
    main()
