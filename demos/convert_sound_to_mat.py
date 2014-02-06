#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import scipy.io
import os

from scikits import audiolab

import cochlea
import mrlib.waves as wv

def convert_sound_to_mat(sound_fname):

    fs = 100e3
    dbspl = 50
    cf = 1000


    ### Read the sound file + resample + scale
    f = audiolab.Sndfile(fname, 'r')
    sound_raw = f.read_frames(f.nframes)
    sound_raw = wv.resample(sound_raw, f.samplerate, fs)
    sound = wv.set_dbspl(s, dbspl)



    ### Run the inner ear model
    anf_trains = cochlea.run_zilany2013(
        sound=sound,
        fs=fs,
        anf_num=(100,75,25),
        cf=cf,
        species='human',
        seed=seed
    )




def main():

    pass


if __name__ == "__main__":
    main()
