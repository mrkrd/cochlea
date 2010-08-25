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


# def bm_init(double fs,
#             np.ndarray[np.float64_t, ndim=1] Ls,
#             np.ndarray[np.float64_t, ndim=1] Rs,
#             np.ndarray[np.float64_t, ndim=1] Ct,
#             np.ndarray[np.float64_t, ndim=1] Rbm,
#             np.ndarray[np.float64_t, ndim=1] Cbm,
#             np.ndarray[np.float64_t, ndim=1] Lbm,
#             double Rh, double Lh):

#     cdef double *Ls_data = <double *>np.PyArray_DATA(Ls)
#     cdef double *Rs_data = <double *>np.PyArray_DATA(Rs)
#     cdef double *Ct_data = <double *>np.PyArray_DATA(Ct)
#     cdef double *Rbm_data = <double *>np.PyArray_DATA(Rbm)
#     cdef double *Cbm_data = <double *>np.PyArray_DATA(Cbm)
#     cdef double *Lbm_data = <double *>np.PyArray_DATA(Lbm)

#     bm_init_c(fs,
#               Ls_data,
#               Rs_data,
#               Ct_data,
#               Rbm_data,
#               Cbm_data,
#               Lbm_data,
#               Rh,
#               Lh)


# def bm_wave(np.ndarray[np.float64_t, ndim=1] signal,
#             np.ndarray[np.float64_t, ndim=1] ampl_corr,
#             np.ndarray[np.float64_t, ndim=1] Abm,
#             np.ndarray[np.float64_t, ndim=1] Cbm):

#     xBM = np.zeros((len(signal), 100))

#     cdef double *signal_data = <double *>np.PyArray_DATA(signal)
#     cdef double *xBM_data = <double *>np.PyArray_DATA(xBM)
#     cdef double *ampl_corr_data = <double *>np.PyArray_DATA(ampl_corr)
#     cdef double *Abm_data = <double *>np.PyArray_DATA(Abm)
#     cdef double *Cbm_data = <double *>np.PyArray_DATA(Cbm)

#     for i in range(len(signal)):
#         bm_wave_c(signal_data[i],
#                   &xBM_data[100*i],
#                   ampl_corr_data,
#                   Abm_data,
#                   Cbm_data)

#     return xBM


import bm_pars

def bm_wave(np.float64_t fs,
            np.ndarray[np.float64_t, ndim=1] signal):


    cdef np.ndarray[np.float64_t] Ls = bm_pars.Ls
    cdef np.ndarray[np.float64_t] Rs = bm_pars.Rs
    cdef np.ndarray[np.float64_t] Ct = bm_pars.Ct
    cdef np.ndarray[np.float64_t] Rbm = bm_pars.Rbm
    cdef np.ndarray[np.float64_t] Cbm = bm_pars.Cbm
    cdef np.ndarray[np.float64_t] Lbm = bm_pars.Lbm
    cdef np.float64_t Rh = bm_pars.Rh
    cdef np.float64_t Lh = bm_pars.Lh
    cdef np.ndarray[np.float64_t] ampl_corr = bm_pars.ampl_corr
    cdef np.ndarray[np.float64_t] Abm = bm_pars.Abm

    cdef np.ndarray[np.float64_t] Z13 = np.zeros(100)
    cdef np.ndarray[np.float64_t] Z32 = np.zeros(100)
    cdef np.ndarray[np.float64_t] Z42 = np.zeros(100)
    cdef np.ndarray[np.float64_t] Z43 = np.zeros(100)
    cdef np.ndarray[np.float64_t] g11 = np.zeros(100)
    cdef np.ndarray[np.float64_t] g12 = np.zeros(100)
    cdef np.ndarray[np.float64_t] g2 = np.zeros(100)
    cdef np.ndarray[np.float64_t] g3 = np.zeros(100)
    cdef np.ndarray[np.float64_t] g41 = np.zeros(100)
    cdef np.ndarray[np.float64_t] g42 = np.zeros(100)
    cdef np.ndarray[np.float64_t] a21 = np.zeros(100)
    cdef np.ndarray[np.float64_t] b14 = np.zeros(100)
    cdef np.ndarray[np.float64_t] b44 = np.zeros(100)
    cdef np.ndarray[np.float64_t] b30 = np.zeros(100)
    cdef np.ndarray[np.float64_t] b33 = np.zeros(100)
    cdef np.ndarray[np.float64_t] b20 = np.zeros(100)
    cdef np.ndarray[np.float64_t] b23 = np.zeros(100)

    cdef Py_ssize_t sec
    cdef Py_ssize_t max_section = 99
    cdef Py_ssize_t i

    cdef np.float64_t R14
    cdef np.float64_t R44
    cdef np.float64_t G33
    cdef np.float64_t G23
    cdef np.float64_t gh
    cdef np.float64_t R_input
    cdef np.float64_t Zhel

    cdef Py_ssize_t signal_len = len(signal)
    cdef np.ndarray[np.float64_t, ndim=2] xbm = np.zeros((signal_len, 100))
    cdef np.float64_t sample
    cdef np.ndarray[np.float64_t] time_slice


    # Init gXX
    for sec in range(max_section,-1,-1):
        if sec == max_section:
            R14 = Rh + (2*fs*Lh)
            gh = Rh / R14

        # Adaptor 4 (series)
        R44 = Rbm[sec] + 2*fs*Lbm[sec] + 1/(2*fs*Cbm[sec])
        g41[sec] = Rbm[sec] / R44
        g42[sec] = 2*fs*Lbm[sec] / R44

        # Adaptor 3 (parallel)
        G33 = 1/R44 + 2*fs*Ct[sec]
        g3[sec] = 1/(G33*R44)

        # Adaptor 2 (parallel)
        G23 = 1/R14 + G33
        g2[sec] = 1 / (R14*G23)

        # Adaptor 1 (series)
        R14 = 1/G23 + Rs[sec] + 2*fs*Ls[sec]
        g11[sec] = 1 / (G23*R14)
        g12[sec] = Rs[sec] / R14

    R_input = R14
    Zhel = 0

    for i in range(signal_len):
        sample = signal[i]
        time_slice = xbm[i]

        # Backward wave
        for sec in range(max_section,-1,-1):
            if sec == max_section:
                a21[sec] = Zhel
            else:
                a21[sec] = -b14[sec+1]

            b44[sec] = -(-Z42[sec] + Z43[sec])

            b30[sec] = -g3[sec]*(Z32[sec] - b44[sec])
            b33[sec]  = Z32[sec] + b30[sec]

            b20[sec] = -g2[sec]*(b33[sec]-a21[sec])
            b23[sec] = b33[sec] + b20[sec]

            b14[sec] = -(b23[sec] - Z13[sec])

        # Forward wave
        for sec in range(100):
            if sec == 0:
                a14 = 2*R_input*sample + b14[0]
            else:
                a14 = -b21

            a10 = a14 - b14[sec]
            b11 = b23[sec] - g11[sec]*a10
            b12 = -g12[sec]*a10
            b13 = -(b11 + b12 + a14)

            b22 = b11 + b20[sec]
            b21 = b22 + b33[sec] - a21[sec]

            b22 = b11 + b20[sec]
            b21 = b22 + b33[sec] - a21[sec]

            b32 = b22 + b30[sec]
            b31 = b32 + Z32[sec] - b44[sec]

            a40 = b31 - b44[sec]
            b41 = -g41[sec]*a40
            b42 = -Z42[sec] - g42[sec]*a40
            b43 = -(b41 + b42 + b31)

            time_slice[sec] = (b43+Z43[sec])*Cbm[sec]/2/Abm[sec]*ampl_corr[sec]

            Z13[sec] = b13
            Z32[sec] = b32
            Z42[sec] = b42
            Z43[sec] = b43


        # Helicotrema
        ah0 = b21 - Zhel
        bh1 = -gh * ah0
        bh2 = -(bh1 + b21)
        Zhel = bh2


    return xbm


def LCR4_init(double fs,
              np.ndarray[np.float64_t, ndim=1] freq_map,
              np.ndarray[np.float64_t, ndim=1] Qmin,
              np.ndarray[np.float64_t, ndim=1] SAT1,
              np.ndarray[np.float64_t, ndim=1] SAT4):

    cdef double *freq_map_data = <double *>np.PyArray_DATA(freq_map)
    cdef double *Qmin_data = <double *>np.PyArray_DATA(Qmin)
    cdef double *SAT1_data = <double *>np.PyArray_DATA(SAT1)
    cdef double *SAT4_data = <double *>np.PyArray_DATA(SAT4)

    LCR4_init_c(fs,
                freq_map_data,
                Qmin_data,
                SAT1_data,
                SAT4_data)



def LCR4(np.ndarray[np.float64_t, ndim=2] xBM,
         np.ndarray[np.float64_t, ndim=1] Qmax,
         np.ndarray[np.float64_t, ndim=1] Qmin):

    xBM = np.array(xBM)         # make a copy, changes are in-place
    cdef double *xBM_data = <double *>np.PyArray_DATA(xBM)
    cdef double *Qmin_data = <double *>np.PyArray_DATA(Qmin)
    cdef double *Qmax_data = <double *>np.PyArray_DATA(Qmax)

    sample_num, section_num = (<object>xBM).shape

    for i in range(sample_num):
        LCR4_c(&xBM_data[i*100],
                Qmax_data,
                Qmin_data,
                0,
                100)

    return xBM



def ihcrp_init(double fs):
    ihcrp_init_c(fs)


def ihcrp(np.ndarray[np.float64_t, ndim=2] xBM,
          np.ndarray[np.float64_t, ndim=1] ciliaGain):

    sample_num, section_num = (<object>xBM).shape

    uIHC = np.zeros_like(xBM)

    cdef double *xBM_data = <double *>np.PyArray_DATA(xBM)
    cdef double *uIHC_data = <double *>np.PyArray_DATA(uIHC)
    cdef double *ciliaGain_data = <double *>np.PyArray_DATA(ciliaGain)

    for i in range(sample_num):
        ihcrp_c(&uIHC_data[section_num*i],
                 &xBM_data[section_num*i],
                 ciliaGain_data)

    return uIHC
