# -*- coding: utf-8 -*-

"""Interoperability with Brian.

"""
from __future__ import division, print_function, absolute_import

import numpy as np

import cochlea
import thorns as th
import thorns.waves as wv

import brian
from brian import mV, second, ms, volt

def main():

    fs = 100e3
    cf = 1e3
    tmax = 50e-3

    sound = wv.ramped_tone(
        fs=fs,
        freq=cf,
        duration=tmax,
        pad=30e-3,
        dbspl=50,
    )

    anf_trains = cochlea.run_zilany2014(
        sound=sound,
        fs=fs,
        cf=cf,
        anf_num=(300, 0, 0),
        seed=0,
        species='cat',
    )

    anfs = cochlea.make_brian_group(anf_trains)

    print(anfs)


    brainstem = make_brainstem_group(100)

    print(brainstem)

    monitor = brian.SpikeMonitor(brainstem)


    net = brian.Network([brainstem, monitor])
    net.run(tmax*second)


    brian.raster_plot(monitor)
    brian.show()



def make_brainstem_group(num):
    """Generate integrate-and-fire neuron group that can be connected to a
    auditory nerve fibers group.

    """
    eqs = '''
    dv/dt = (ge-(v+49*mV))/(2*ms) : volt
    dge/dt = -ge/(0.5*ms) : volt
    '''

    group = brian.NeuronGroup(num, eqs, threshold=-50*mV, reset=-60*mV)
    group.v = -60*mV

    return group



if __name__ == "__main__":
    main()
