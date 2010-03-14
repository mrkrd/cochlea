import numpy as np
cimport numpy as np
from stdlib cimport malloc
import scipy.signal as dsp

cdef extern from "stdlib.h":
    void *memcpy(void *str1, void *str2, size_t n)

np.import_array()

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
    signal_arr = np.PyArray_SimpleNewFromData(1, [k],
                                              np.NPY_DOUBLE,
                                              <void*>signal)

    # Filter + downsample
    b = dsp.firwin(q*2, 1./q, window='hamming')
    filtered = dsp.lfilter(b, 1., signal_arr)
    decimated = np.array(filtered[::q])


    # Copy data to output array
    cdef double *decimated_ptr = <double *>np.PyArray_DATA(decimated)
    cdef double *out_ptr = <double *>malloc(len(decimated)*sizeof(double))
    memcpy(out_ptr, decimated_ptr, len(decimated)*sizeof(double))

    return out_ptr
