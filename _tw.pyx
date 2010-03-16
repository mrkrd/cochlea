import numpy as np
cimport numpy as np


# cdef extern from "numpy/arrayobject.h":
#     void import_array()



cdef extern from "bm_wave.h":
    double bm_init_c(double f_s,
                     double *Ls,
                     double *Rs,
                     double *Ct,
                     double *Rbm,
                     double *Cbm,
                     double *Lbm,
                     double Rh,
                     double Lh)
    void bm_wave_c(double input,
                   double *x_BM,
                   double *ampl_corr,
                   double *Abm,
                   double *Cbm)


cdef extern from "LCR4.h":
    void LCR4_init_c(double f_s,
                     double *freq_map,
                     double *Qmin,
                     double *SAT1,
                     double *SAT4)

    void LCR4_c(double *xBM,
                double *Qmax,
                double *Qmin,
                int first_sec,
                int last_sec)


cdef extern from "ihcrp.h":
    double ihcrp_init_c(double f_s)
    void ihcrp_c(double *uIHC,
                 double *xBM,
                 double *ciliaCouplingGain)



# import_array()


# TODO: type check
def bm_init(double fs, Ls, Rs, Ct, Rbm, Cbm, Lbm, double Rh, double Lh):

    cdef double *Ls_data = <double *>np.PyArray_DATA(Ls)
    cdef double *Rs_data = <double *>np.PyArray_DATA(Rs)
    cdef double *Ct_data = <double *>np.PyArray_DATA(Ct)
    cdef double *Rbm_data = <double *>np.PyArray_DATA(Rbm)
    cdef double *Cbm_data = <double *>np.PyArray_DATA(Cbm)
    cdef double *Lbm_data = <double *>np.PyArray_DATA(Lbm)

    bm_init_c(fs,
              Ls_data,
              Rs_data,
              Ct_data,
              Rbm_data,
              Cbm_data,
              Lbm_data,
              Rh,
              Lh)


def bm_wave(signal, ampl_corr, Abm, Cbm):

    xBM = np.zeros((len(signal), 100))

    cdef double *signal_data = <double *>np.PyArray_DATA(signal)
    cdef double *xBM_data = <double *>np.PyArray_DATA(xBM)
    cdef double *ampl_corr_data = <double *>np.PyArray_DATA(ampl_corr)
    cdef double *Abm_data = <double *>np.PyArray_DATA(Abm)
    cdef double *Cbm_data = <double *>np.PyArray_DATA(Cbm)

    for i in range(len(signal)):
        bm_wave_c(signal_data[i],
                  &xBM_data[100*i],
                  ampl_corr_data,
                  Abm_data,
                  Cbm_data)

    return xBM



def LCR4_init(double fs, freq_map, Qmin, SAT1, SAT4):

    cdef double *freq_map_data = <double *>np.PyArray_DATA(freq_map)
    cdef double *Qmin_data = <double *>np.PyArray_DATA(Qmin)
    cdef double *SAT1_data = <double *>np.PyArray_DATA(SAT1)
    cdef double *SAT4_data = <double *>np.PyArray_DATA(SAT4)

    LCR4_init_c(fs,
                freq_map_data,
                Qmin_data,
                SAT1_data,
                SAT4_data)



def LCR4(xBM, Qmax, Qmin):

    xBM = np.array(xBM)         # don't do inplace calculations
    cdef double *xBM_data = <double *>np.PyArray_DATA(xBM)
    cdef double *Qmin_data = <double *>np.PyArray_DATA(Qmin)
    cdef double *Qmax_data = <double *>np.PyArray_DATA(Qmax)

    sample_num, section_num = xBM.shape

    for i in range(sample_num):
        LCR4_c(&xBM_data[i*100],
                Qmin_data,
                Qmax_data,
                0,
                100)

    return xBM



def ihcrp_init(double fs):
    ihcrp_init_c(fs)


def ihcrp(xBM, ciliaGain):

    sample_num, section_num = xBM.shape

    uIHC = np.zeros_like(xBM)

    cdef double *xBM_data = <double *>np.PyArray_DATA(xBM)
    cdef double *uIHC_data = <double *>np.PyArray_DATA(uIHC)
    cdef double *ciliaGain_data = <double *>np.PyArray_DATA(ciliaGain)

    for i in range(sample_num):
        ihcrp_c(&uIHC_data[section_num*i],
                 &xBM_data[section_num*i],
                 ciliaGain_data)

    return uIHC
