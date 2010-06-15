import numpy as np
from stdlib cimport malloc
import scipy.signal as dsp
import ffGn_module

cimport numpy as np


cdef extern from "stdlib.h":
    void *memcpy(void *str1, void *str2, size_t n)


cdef extern from "catmodel.h":
    void IHCAN(double *px,
               double cf,
               int nrep,
               double tdres,
               int totalstim,
               double cohc,
               double cihc,
               double *ihcout)
    double Synapse(double *ihcout,
                   double tdres,
                   double cf,
                   int totalstim,
                   int nrep,
                   double spont,
                   double implnt,
                   double sampFreq,
                   double *synouttmp,
                   int with_ffGn)
    int SpikeGenerator(double *synouttmp,
                       double tdres,
                       int totalstim,
                       int nrep,
                       double *sptime)


cdef extern from "Python.h":
    ctypedef int Py_intptr_t


cdef extern from "numpy/arrayobject.h":
    ctypedef Py_intptr_t npy_intp
    object PyArray_SimpleNewFromData(int nd, npy_intp* dims, int typenum, void* data)
    void import_array()


import_array()




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
    # b = dsp.firwin(2*q, 1./q, window='hamming')
    # filtered = dsp.lfilter(b, 1., signal_arr)
    # decimated = np.array(filtered[::q])

    new_len = len(signal_arr) // q
    signal_arr = signal_arr[0:new_len*q]
    decimated = signal_arr.reshape((new_len, q))
    b = dsp.firwin(q, 1./q, window='hamming')
    decimated = np.sum(decimated*b, axis=1)


    # Copy data to output array
    cdef double *decimated_ptr = <double *>np.PyArray_DATA(decimated)
    cdef double *out_ptr = <double *>malloc(len(decimated)*sizeof(double))
    memcpy(out_ptr, decimated_ptr, len(decimated)*sizeof(double))

    return out_ptr


cdef public double* ffGn(int N, double tdres, double Hinput, double mu, bint with_ffGn):
    """ ffGn.py wrapper """

    if with_ffGn:
        a = ffGn_module.ffGn(N, tdres, Hinput, mu)
    else:
        a = np.zeros(N)

    # Copy data to output array
    cdef double *ptr = <double *>np.PyArray_DATA(a)
    cdef double *out_ptr = <double *>malloc(len(a)*sizeof(double))
    memcpy(out_ptr, ptr, len(a)*sizeof(double))

    return out_ptr



# TODO: type check
def run_ihc(np.ndarray[np.float64_t, ndim=1] signal,
            double cf,
            double fs,
            int cohc=1,
            int cihc=1):
    """
    Run BM / IHC model.

    signal: input sound in uPa
    cf: characteristic frequency
    fs: sampling frequency
    cohc, cihc: degeneration parameters for IHC and OHC cells

    return: IHC receptor potential

    """
    assert (cf > 80) and (cf < 40e3), "Wrong CF: 80 < cf < 40e3"
    assert (fs >= 100e3) and (fs <= 500e3), "Wrong Fs: 100e3 <= fs <= 500e3"
    assert (cohc >= 0) and (cohc <= 1), "0 <= cohc <= 1"
    assert (cihc >= 0) and (cihc <= 1), "0 <= cihc <= 1"

    # uPa -> Pa
    # Compatibility with DSAM
    signal = signal * 1e-6


    # Input sound
    cdef double *signal_data = <double *>np.PyArray_DATA(signal)

    # Output IHC voltage
    ihcout = np.zeros(len(signal))
    cdef double *ihcout_data = <double *>np.PyArray_DATA(ihcout)


    IHCAN(signal_data, cf, 1, 1.0/fs, len(signal), cohc, cihc, ihcout_data);


    return ihcout




def run_synapse(np.ndarray[np.float64_t, ndim=1] vihc,
                double cf,
                double fs,
                anf_type='hsr',
                powerlaw_implnt='actual',
                with_ffGn=True):
    """
    Run synapse simulation.

    vihc: IHC receptor potential
    cf: characteristic frequency
    anf_type: auditory nerve fiber type ('hsr', 'msr' or 'lsr')
    powerlaw_implnt: implementation of the powerlaw ('actual', 'approx')
    with_ffGn: enable/disable factorial Gauss noise generator

    return: PSTH from ANF
    """
    assert (cf > 80) and (cf < 40e3), "Wrong CF: 80 < cf < 40e3"
    assert (fs >= 100e3) and (fs <= 500e3), "Wrong Fs: 100e3 <= fs <= 500e3"
    assert anf_type in ['hsr', 'msr', 'lsr'], "anf_type not hsr/msr/lsr"
    assert powerlaw_implnt in ['actual', 'approx'], "powerlaw_implnt not actual/approx"

    spont = {'hsr': 100.,
             'msr': 5.,
             'lsr': 0.1}

    implnt_map = {'actual': 1,
                  'approx': 0}


    # Input IHC voltage
    cdef double *vihc_data = <double *>np.PyArray_DATA(vihc)


    # Output synapse data (spiking probabilities)
    synout = np.zeros_like(vihc)
    cdef double *synout_data = <double *>np.PyArray_DATA(synout)


    # Run synapse model
    Synapse(vihc_data,                   # px
            1.0/fs,                      # tdres
            cf,                          # cf
            len(vihc),                   # totalstim
            1,                           # nrep
            spont[anf_type],             # spont
            implnt_map[powerlaw_implnt], # implnt
            10e3,                        # sampFreq
            synout_data,                 # synouttmp
            with_ffGn)                   # with_ffGn

    return synout





def run_spike_generator(np.ndarray[np.float64_t, ndim=1] synout,
                        double fs):
    """
    Run spike generator.

    synout: synapse output (sp/s)
    fs: sampling frequency

    return: sptime

    """
    # Input IHC voltage
    cdef double *synout_data = <double *>np.PyArray_DATA(synout)

    # Output spikes (signal)
    sptime = np.zeros(np.ceil(len(synout)/0.00075/fs))
    cdef double *sptime_data = <double *>np.PyArray_DATA(sptime)


    # Run synapse model
    SpikeGenerator(synout_data,  # synouttmp
                   1./fs,        # tdres
                   len(synout),  # totalstim
                   1,            # nprep
                   sptime_data)  # sptime

    return sptime





def set_dbspl(dB, signal):
    p0 = 2e-5                   # Pa
    squared = signal**2
    rms = np.sqrt( np.sum(squared) / len(signal) )

    if rms == 0:
        r = 0
    else:
        r = 10**(dB / 20.0) * p0 / rms

    return signal * r * 1e6     # uPa

