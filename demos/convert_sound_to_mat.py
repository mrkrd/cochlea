#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import scipy.io
import os
import sys
import logging

from scikits import audiolab

import cochlea
import mrlib.waves as wv


def convert_sound_to_mat(sound_fname):

    print("Processing " + sound_fname)

    fs = 100e3
    dbspl = 50
    cf = 1000


    ### Read the sound file + resample + scale
    f = audiolab.Sndfile(sound_fname, 'r')
    sound_raw = f.read_frames(f.nframes)
    sound_raw = wv.resample(sound_raw, f.samplerate, fs)
    sound = wv.set_dbspl(sound_raw, dbspl)



    ### Run the inner ear model
    anf_trains = cochlea.run_zilany2013(
        sound=sound,
        fs=fs,
        anf_num=(100,75,25),
        cf=cf,
        species='human',
        seed=0
    )



    ### Save spikes to matlab file
    trains = anf_trains.to_records()
    mat_fname = os.path.splitext(sound_fname)[0]
    mdict = {'trains': trains}

    scipy.io.savemat(
        mat_fname,
        mdict,
        do_compression=True
    )



def main():

    if sys.argv < 2:
        print("Usage: python convert_sounds_to_mat.py a.wav b.wav c.sph")
        exit()

    map(convert_sound_to_mat, sys.argv[1:])

    # TODO: add some sounds to the repository


if __name__ == "__main__":
    main()
