import numpy as np
from stdlib cimport malloc
import scipy.signal as dsp
import ffGn_module

cimport numpy as np


cdef extern from "stdlib.h":
    void *memcpy(void *str1, void *str2, size_t n)


cdef extern from "catmodel.h":
    void IHCAN(double *px, double cf, int nrep, double tdres, int totalstim,
               double cohc, double cihc, double *ihcout)
    void SingleAN(double *px, double cf, int nrep, double tdres, int totalstim,
                  double fibertype, double implnt, double *synout, double *psth)


cdef extern from "Python.h":
    ctypedef int Py_intptr_t


cdef extern from "numpy/arrayobject.h":
    ctypedef Py_intptr_t npy_intp
    object PyArray_SimpleNewFromData(int nd, npy_intp* dims, int typenum, void* data)
    void import_array()


import_array()



cdef public int is_pycat_initialized = 0


cdef public double* generate_random_numbers(long length):
    arr = np.random.rand(length)

    cdef double *data_ptr = <double *>np.PyArray_DATA(arr)

    cdef double *out_ptr = <double *>malloc(length * sizeof(double))
    memcpy(out_ptr, data_ptr, length*sizeof(double))

    return out_ptr




cdef public double* decimate(int k, double *signal, int q):
    """ Decimate signal.

    k: lenght of signal
    signal: pointer to the signal
    q: decimation factor

    """
    # signal_arr will not own the data, signal's array has to be freed
    # after return from this function
    signal_arr = PyArray_SimpleNewFromData(1, [k],
                                           np.NPY_DOUBLE,
                                           <void *>signal)

    # Filter + downsample
    b = dsp.firwin(2*q, 1./q, window='hamming')
    filtered = dsp.lfilter(b, 1., signal_arr)
    decimated = np.array(filtered[::q])


    # Copy data to output array
    cdef double *decimated_ptr = <double *>np.PyArray_DATA(decimated)
    cdef double *out_ptr = <double *>malloc(len(decimated)*sizeof(double))
    memcpy(out_ptr, decimated_ptr, len(decimated)*sizeof(double))

    return out_ptr


cdef public double* ffGn(int N, double tdres, double Hinput, double mu):
    """ ffGn.py wrapper """

    a = ffGn_module.ffGn(N, tdres, Hinput, mu)

    # Copy data to output array
    cdef double *ptr = <double *>np.PyArray_DATA(a)
    cdef double *out_ptr = <double *>malloc(len(a)*sizeof(double))
    memcpy(out_ptr, ptr, len(a)*sizeof(double))

    return out_ptr



# TODO: type check
def run_ihc(signal, cf, fs, cohc=1, cihc=1, verbose=False):
    """
    signal: input sound in uPa
    cf: characteristic frequency
    fs: sampling frequency
    cohc, cihc: degeneration parameters for IHC and OHC cells

    return: IHC receptor potential
    """
    assert isinstance(signal, np.ndarray) and (signal.ndim == 1)
    assert (cf > 80) and (cf < 40e3)
    assert (fs >= 100e3) and (fs <= 500e3)
    assert (cohc >= 0) and (cohc <= 1)
    assert (cihc >= 0) and (cihc <= 1)

    # uPa -> Pa
    # Compatibility with DSAM
    signal = signal * 1e-6

    if verbose:
        print "IHC@", cf


    # Input sound
    cdef double *signal_data = <double *>np.PyArray_DATA(signal)

    # Output IHC voltage
    ihcout = np.zeros(len(signal))
    cdef double *ihcout_data = <double *>np.PyArray_DATA(ihcout)


    IHCAN(signal_data, cf, 1, 1.0/fs, len(signal), cohc, cihc, ihcout_data);


    return ihcout




def run_synapse(vihc, cf, fs, anf_type='hsr', powerlaw_implnt='actual',
                verbose=False):
    """
    vihc: IHC receptor potential
    cf: characteristic frequency
    anf_type: auditory nerve fiber type ('hsr', 'msr' or 'lsr')
    powerlaw_implnt: implementation of the powerlaw ('actual', 'approx')

    return: PSTH from ANF
    """
    assert isinstance(vihc, np.ndarray) and (vihc.ndim == 1)
    assert (cf > 80) and (cf < 40e3)
    assert (fs >= 100e3) and (fs <= 500e3)
    assert anf_type in ['hsr', 'msr', 'lsr']
    assert powerlaw_implnt in ['actual', 'approx']

    anf_map = {'hsr': 3,
               'msr': 2,
               'lsr': 1}

    implnt_map = {'actual': 1,
                  'approx': 0}

    if verbose:
        print "Syn@", cf

    # Input IHC voltage
    cdef double *vihc_data = <double *>np.PyArray_DATA(vihc)

    # Output spikes (signal)
    spikes = np.zeros_like(vihc)
    cdef double *spikes_data = <double *>np.PyArray_DATA(spikes)

    # Output synapse data (spiking probabilities)
    synout = np.zeros_like(vihc)
    cdef double *synout_data = <double *>np.PyArray_DATA(synout)


    # Run synapse model
    SingleAN(vihc_data, cf, 1, 1.0/fs, len(vihc),
             anf_map[anf_type], implnt_map[powerlaw_implnt],
             synout_data, spikes_data);

    return spikes



def set_dbspl(dB, signal):
    p0 = 2e-5                   # Pa
    squared = signal**2
    rms = np.sqrt( np.sum(squared) / len(signal) )

    if rms == 0:
        r = 0
    else:
        r = 10**(dB / 20.0) * p0 / rms

    return signal * r * 1e6     # uPa

