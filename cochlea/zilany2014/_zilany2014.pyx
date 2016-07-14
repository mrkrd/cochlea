# Copyright 2013-2014 Marek Rudnicki

# This file is part of cochlea.

# cochlea is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# cochlea is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with cochlea.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division, print_function, absolute_import

import numpy as np
from libc.stdlib cimport malloc
from . import util
import scipy.signal as dsp

cimport numpy as np

cdef extern from "stdlib.h":
    void *memcpy(void *str1, void *str2, size_t n)


cdef extern from "model_IHC.h":
    void IHCAN(
        double *px,
        double cf,
        int nrep,
        double tdres,
        int totalstim,
        double cohc,
        double cihc,
        int species,
        double *ihcout
    )

cdef extern from "model_Synapse.h":
    double Synapse(
        double *ihcout,
        double tdres,
        double cf,
        int totalstim,
        int nrep,
        double spont,
        double noiseType,
        double implnt,
        double sampFreq,
        double *synouttmp
    )
    int SpikeGenerator(
        double *synouttmp,
        double tdres,
        int totalstim,
        int nrep,
        double *sptime
    )


cdef extern from "Python.h":
    ctypedef int Py_intptr_t


cdef extern from "numpy/arrayobject.h":
    ctypedef Py_intptr_t npy_intp
    object PyArray_SimpleNewFromData(
        int nd,
        npy_intp* dims,
        int typenum,
        void* data
    )



np.import_array()


def run_ihc(
        np.ndarray[np.float64_t, ndim=1] signal,
        double cf,
        double fs,
        species='cat',
        double cohc=1.,
        double cihc=1.
):
    """Run middle ear filter, BM filters and IHC model.

    Parameters
    ----------
    signal : array_like
        Output of the middle ear filter in Pascal.
    cf : float
        Characteristic frequency in Hz.
    fs : float
        Sampling frequency in Hz.
    species : {'cat', 'human', 'human_glasberg1990'}
        Species.
    cihc, cohc : float
        Degeneration parameters for IHC and OHC cells.

    Returns
    -------
    array_like
        IHC receptor potential.

    """
    if species == 'cat':
        assert (cf > 124.9) and (cf < 40e3), "Wrong CF: 125 <= cf < 40e3, CF = %s"%str(cf)
    elif 'human' in species:
        assert (cf > 124.9) and (cf < 20001.), "Wrong CF: 125 <= cf <= 20e3, CF = %s"%str(cf)

    assert (fs >= 100e3) and (fs <= 500e3), "Wrong Fs: 100e3 <= fs <= 500e3"
    assert (cohc >= 0) and (cohc <= 1), "0 <= cohc <= 1"
    assert (cihc >= 0) and (cihc <= 1), "0 <= cihc <= 1"


    species_map = {
        'cat': 1,
        'human': 2,
        'human_glasberg1990': 3,
    }

    # Input sound
    if not signal.flags['C_CONTIGUOUS']:
        signal = signal.copy(order='C')
    cdef double *signal_data = <double *>np.PyArray_DATA(signal)

    # Output IHC voltage
    ihcout = np.zeros( len(signal) )
    cdef double *ihcout_data = <double *>np.PyArray_DATA(ihcout)


    IHCAN(
        signal_data,
        cf,
        1,
        1.0/fs,
        len(signal),
        cohc,
        cihc,
        species_map[species],
        ihcout_data
    )


    return ihcout




def run_synapse(
        np.ndarray[np.float64_t, ndim=1] vihc,
        double fs,
        double cf,
        anf_type='hsr',
        powerlaw='actual',
        ffGn=True
):
    """Run synapse simulation.

    vihc: IHC receptor potential
    cf: characteristic frequency
    anf_type: auditory nerve fiber type ('hsr', 'msr' or 'lsr')
    powerlaw: implementation of the powerlaw ('actual', 'approximate')
    ffGn: enable/disable factorial Gauss noise generator

    return: PSTH from ANF

    """
    assert (cf > 79.9) and (cf < 40e3), "Wrong CF: 80 <= cf < 40e3, CF = %s"%str(cf)
    assert (fs >= 100e3) and (fs <= 500e3), "Wrong Fs: 100e3 <= fs <= 500e3"
    assert anf_type in ['hsr', 'msr', 'lsr'], "anf_type not hsr/msr/lsr"
    assert powerlaw in ['actual', 'approximate'], "powerlaw not actual/approximate"

    spont = {
        'lsr': 0.1,
        'msr': 4.0,
        'hsr': 100.0,
    }

    powerlaw_map = {
        'actual': 1,
        'approximate': 0
    }

    if ffGn:
        noise_type = 1.
    else:
        noise_type = 0.


    # Input IHC voltage
    if not vihc.flags['C_CONTIGUOUS']:
        vihc = vihc.copy(order='C')
    cdef double *vihc_data = <double *>np.PyArray_DATA(vihc)


    # Output synapse data (spiking probabilities)
    synout = np.zeros_like(vihc)
    cdef double *synout_data = <double *>np.PyArray_DATA(synout)


    # Run synapse model
    Synapse(
        vihc_data,                   # ihcout
        1.0/fs,                      # tdres
        cf,                          # cf
        len(vihc),                   # totalstim
        1,                           # nrep
        spont[anf_type],             # spont
        noise_type,                  # noiseType
        powerlaw_map[powerlaw],      # implnt
        10e3,                        # sampFreq
        synout_data                  # synouttmp
    )

    return synout


def run_spike_generator(
        np.ndarray[np.float64_t, ndim=1] synout,
        double fs
):
    """Run spike generator.

    synout: synapse output
    fs: sampling frequency

    return: sptime

    """
    # Input IHC voltage
    if not synout.flags['C_CONTIGUOUS']:
        synout = synout.copy(order='C')
    cdef double *synout_data = <double *>np.PyArray_DATA(synout)

    # Output spikes (signal)
    sptimes = np.zeros(int(np.ceil(len(synout)/0.00075/fs)))
    cdef double *sptimes_data = <double *>np.PyArray_DATA(sptimes)


    # Run synapse model
    SpikeGenerator(
        synout_data,            # synouttmp
        1./fs,                  # tdres
        len(synout),            # totalstim
        1,                      # nprep
        sptimes_data            # sptime
    )

    spikes = np.array(sptimes[sptimes != 0])

    return spikes



cdef public double* generate_random_numbers(long length):
    arr = np.random.rand(length)

    if not arr.flags['C_CONTIGUOUS']:
        arr = arr.copy(order='C')

    cdef double *data_ptr = <double *>np.PyArray_DATA(arr)
    cdef double *out_ptr = <double *>malloc(length * sizeof(double))
    memcpy(out_ptr, data_ptr, length*sizeof(double))

    return out_ptr




cdef public double* decimate(
    int k,
    double *signal,
    int q
):
    """Decimate a signal

    k: number of samples in signal
    signal: pointer to the signal
    q: decimation factor

    This implementation was inspired by scipy.signal.decimate.

    """
    # signal_arr will not own the data, signal's array has to be freed
    # after return from this function
    signal_arr = PyArray_SimpleNewFromData(
        1,                      # nd
        [k],                    # dims
        np.NPY_DOUBLE,          # typenum
        <void *>signal          # data
    )


    # resampled = dsp.resample(
    #     signal_arr,
    #     len(signal_arr) // q
    # )


    b = dsp.firwin(q+1, 1./q, window='hamming')
    a = [1.]

    filtered = dsp.filtfilt(
        b=b,
        a=a,
        x=signal_arr
    )

    resampled = filtered[::q]


    if not resampled.flags['C_CONTIGUOUS']:
        resampled = resampled.copy(order='C')


    # Copy data to output array
    cdef double *resampled_ptr = <double *>np.PyArray_DATA(resampled)
    cdef double *out_ptr = <double *>malloc(len(resampled)*sizeof(double))
    memcpy(out_ptr, resampled_ptr, len(resampled)*sizeof(double))

    return out_ptr


cdef public double* ffGn(int N, double tdres, double Hinput, double noiseType, double mu):
    """util.ffGn() wrapper"""

    a = util.ffGn(N, tdres, Hinput, noiseType, mu)

    if not a.flags['C_CONTIGUOUS']:
        a = a.copy(order='C')

    # Copy data to output array
    cdef double *ptr = <double *>np.PyArray_DATA(a)
    cdef double *out_ptr = <double *>malloc(len(a)*sizeof(double))
    memcpy(out_ptr, ptr, len(a)*sizeof(double))

    return out_ptr
