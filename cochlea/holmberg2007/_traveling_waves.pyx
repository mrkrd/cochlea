from __future__ import division

import numpy as np
cimport numpy as np


cdef extern from "math.h":
    double sqrt(double x)
    double exp(double x)
    double fabs(double x)





import bm_pars

def run_bm_wave(
        np.ndarray[np.float64_t, ndim=1] signal,
        np.float64_t fs,
):

    assert fs == 48e3

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

    cdef np.float64_t gh, R_input, Zhel
    cdef np.float64_t R14, R44, G33, G23
    cdef np.float64_t a10, b11, b12, b13, a40, b41, b42, b43
    cdef np.float64_t b21, b22, b31, b32, ah0, bh1, bh2, a14


    cdef Py_ssize_t sec
    cdef Py_ssize_t i

    cdef np.ndarray[np.float64_t, ndim=2] xbm = np.zeros((len(signal), 100))
    cdef np.float64_t sample
    cdef np.ndarray[np.float64_t] time_slice


    cdef Py_ssize_t sections = 100


    # Init gXX
    for sec in range(sections-1,-1,-1):
        if sec == sections-1:
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

    # Iterate over time
    for i in range(len(signal)):
        sample = signal[i]
        time_slice = xbm[i]

        # Backward wave
        for sec in range(sections-1,-1,-1):
            if sec == sections-1:
                a21[sec] = Zhel
            else:
                a21[sec] = -b14[sec+1]

            b44[sec] = -(-Z42[sec] + Z43[sec])

            b30[sec] = -g3[sec]*(Z32[sec] - b44[sec])
            b33[sec] = Z32[sec] + b30[sec]

            b20[sec] = -g2[sec]*(b33[sec]-a21[sec])
            b23[sec] = b33[sec] + b20[sec]

            b14[sec] = -(b23[sec] - Z13[sec])

        # Forward wave
        for sec in range(sections):
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


    flipped = np.fliplr(xbm)
    return flipped



def run_lcr4(xbm, fs, sections=None):

    assert fs == 48e3

    if sections is None:
        sections = range(100)
    elif isinstance(sections, int):
        sections = [sections]

    if xbm.ndim == 1:
        xbm = xbm.reshape(xbm.shape+(1,))

    assert xbm.shape[1] == len(sections)

    xbm_out = []
    for i in range(len(sections)):
        sec = sections[i]
        xbm_slice = xbm[:,i]
        xbm_out.append( _run_single_lcr4(fs, xbm_slice, sec) )

    return np.array(xbm_out).T

def _run_single_lcr4(np.float64_t fs,
                     np.ndarray[np.float64_t] xbm,
                     Py_ssize_t sec):

    assert fs == 48000

    sec = 99 - sec         # DSAM convention (low channel -> low freq)

    cdef np.ndarray[np.float64_t] freq_map = bm_pars.freq_map
    cdef np.ndarray[np.float64_t] Qmin = bm_pars.Qmin
    cdef np.ndarray[np.float64_t] Qmax = bm_pars.Qmax
    cdef np.ndarray[np.float64_t] SAT1 = bm_pars.SAT1
    cdef np.ndarray[np.float64_t] SAT4 = bm_pars.SAT4

    cdef np.float64_t T4C, T4L, T5C, T5L, T6C, T6L
    cdef np.float64_t T7C, T7L, Rtest, Rc_res, Rl_res
    cdef np.float64_t Z_CLP1,Z_CLP2,Z_CLP3,Z_CLP4,g_CLP
    cdef np.float64_t p_SAT1,p_SAT2,p_SAT3,p_SAT4
    cdef np.float64_t Qd0, Qd1, Qd2, Qd3

    cdef np.float64_t R1, R2, tau

    cdef np.float64_t V2, V3, V4, V5
    cdef np.float64_t X2, X3, X4, X5
    cdef np.float64_t Rr11, Rr13, Gr21, Gr23, gr1, gr2
    cdef np.float64_t b11, b12, b13, b20, b21, b22, b23, a0

    cdef np.float64_t rect, rect_LP, b2

    cdef Py_ssize_t i

    tau = 1 / (2*np.pi*0.8e3)
    R2 = 1 / (2*fs*tau)
    R1 = 1

    g_CLP = (R2-R1) / (R2+R1)

    Rc_res = Qmin[sec] / (2*fs)
    Rl_res = Rc_res * (fs/(np.pi*freq_map[sec]))**2

    T4C = 0
    T4L = 0
    T5C = 0
    T5L = 0
    T6C = 0
    T6L = 0
    T7C = 0
    T7L = 0

    Qd0 = Qmin[sec]
    Qd1 = Qmin[sec]
    Qd2 = Qmin[sec]
    Qd3 = Qmin[sec]

    p_SAT1 = SAT1[sec]
    p_SAT4 = SAT4[sec]

    p_SAT2 = 10**(np.log10(p_SAT4) + np.log10(p_SAT1/p_SAT4)*2/3)
    p_SAT3 = 10**(np.log10(p_SAT4) + np.log10(p_SAT1/p_SAT4)/3)

    Z_CLP1 = 0
    Z_CLP2 = 0
    Z_CLP3 = 0
    Z_CLP4 = 0

    xbm_out = np.zeros_like(xbm)

    for i in range(len(xbm)):

        ### Resonator 1
        Rr11 = sqrt(Rl_res*Rc_res) / Qd0
        Rr13 = Rr11 + Rl_res

        Gr21 = 1 / Rr13
        Gr23 = Gr21 + 1 / Rc_res

        gr1 = Rr11 / Rr13
        gr2 = Gr21 / Gr23

        ### Filtering
        b13 = -(xbm[i] + T4L)
        b20 = -gr2*(T4C - b13)
        b23 = b20 + T4C

        b22 = b20+b23
        b21 = b22+T4C-b13
        a0 = b21-b13
        b11 = xbm[i]-gr1*a0
        b12 = -(b11+b21)

        X2 = (b22+T4C)/2
        V2 = (T4C-b22)/(2*Rc_res)

        T4L = -b12
        T4C = b22

        rect = 1- (2. / (1+ exp( -p_SAT1*fabs(X2)) ) -1)
        b2 = rect*(1 + g_CLP) - Z_CLP1*g_CLP
        rect_LP = (b2 + Z_CLP1)/2.
        Z_CLP1 = b2

        Qd0 = (Qmax[sec]-Qmin[sec])*rect_LP+Qmin[sec]

        ### Resonator 2
        # Calculating port values from old Qd-value
        Rr11 = sqrt(Rl_res*Rc_res) / Qd1 # sqrt(L/C)/Q
        Rr13 = Rr11 + Rl_res

        Gr21 = 1/Rr13
        Gr23 = Gr21 + 1/Rc_res

        gr1 = Rr11/Rr13
        gr2 = Gr21/Gr23

        ### Filtering
        b13 = -(X2 + T5L)
        b20 = -gr2*(T5C - b13)
        b23 = b20 + T5C

        b22 = b20+b23
        b21 = b22+T5C-b13
        a0 = b21-b13
        b11 = X2-gr1*a0
        b12 = -(b11+b21)

        X3 = (b22+T5C)/2.
        V3 = (T5C-b22)/(2.*Rc_res)

        T5L = -b12
        T5C = b22

        rect = 1- (2. / (1+exp( -p_SAT2*fabs(X3)) ) -1)
        b2 = rect*(1 + g_CLP) - Z_CLP2*g_CLP
        rect_LP = (b2 + Z_CLP2)/2.
        Z_CLP2 = b2

        Qd1 = (Qmax[sec]-Qmin[sec])*rect_LP+Qmin[sec]


        ### Resonator 3
        Rr11 = sqrt(Rl_res*Rc_res)/Qd2
        Rr13 = Rr11 + Rl_res

        Gr21 = 1/Rr13
        Gr23 = Gr21 + 1/Rc_res

        gr1 = Rr11/Rr13
        gr2 = Gr21/Gr23

        ### Filtering
        b13 = -(X3+T6L)
        b20 = -gr2*(T6C-b13)
        b23 = b20+T6C

        b22 = b20+b23
        b21 = b22+T6C-b13
        a0 = b21-b13
        b11 = X3-gr1*a0
        b12=-(b11+b21)

        X4 = (b22+T6C)/2.
        V4 = (T6C-b22)/(2.*Rc_res)


        T6L = -b12
        T6C = b22

        rect = 1- (2. / (1+exp( -p_SAT3*fabs(X4)) ) -1)
        b2 = rect*(1 + g_CLP) - Z_CLP3*g_CLP
        rect_LP = (b2 + Z_CLP3)/2.
        Z_CLP3 = b2

        Qd2 = (Qmax[sec]-Qmin[sec])*rect_LP+Qmin[sec]

        ### Resonator 4
        Rr11 = sqrt(Rl_res*Rc_res)/Qd3
        Rr13 = Rr11 + Rl_res

        Gr21 = 1/Rr13
        Gr23 = Gr21 + 1/Rc_res

        gr1 = Rr11/Rr13
        gr2 = Gr21/Gr23

        ### Filtering
        b13 = -(X4+T7L)
        b20 = -gr2*(T7C-b13)
        b23 = b20+T7C

        b22 = b20+b23
        b21 = b22+T7C-b13
        a0 = b21-b13
        b11 = X4-gr1*a0
        b12 = -(b11+b21)

        X5 = (b22+T7C)/2.
        V5 = (T7C-b22)/(2.*Rc_res)

        T7L = -b12
        T7C = b22

        rect = 1 - (2 / (1 + exp( -p_SAT4*fabs(X5)) ) - 1)
        b2 = rect*(1 + g_CLP) - Z_CLP4*g_CLP
        rect_LP = (b2 + Z_CLP4) / 2
        Z_CLP4 = b2

        Qd3 = (Qmax[sec]-Qmin[sec])*rect_LP+Qmin[sec]



        xbm_out[i] = X5

    return xbm_out



def run_ihcrp(xbm, fs, sections=None):

    assert fs == 48e3

    if sections is None:
        sections = range(100)
    elif isinstance(sections, int):
        sections = [sections]

    if xbm.ndim == 1:
        xbm = xbm.reshape(xbm.shape+(1,))

    assert xbm.shape[1] == len(sections)

    uihc = []
    for i in range(len(sections)):
        sec = sections[i]
        xbm_slice = xbm[:,i]
        uihc.append( _run_single_ihcrp(fs, xbm_slice, sec) )

    return np.array(uihc).T



def _run_single_ihcrp(np.float64_t fs,
                      np.ndarray[np.float64_t] xbm,
                      Py_ssize_t sec):

    assert fs == 48000

    sec = 99 - sec         # DSAM convention (low channel -> low freq)

    cdef np.ndarray[np.float64_t] uIHC = np.zeros_like(xbm)
    cdef np.ndarray[np.float64_t] ciliaCouplingGain = bm_pars.ciliaGain

    cdef np.float64_t p_endocochlearPot_Et = 0.1
    cdef np.float64_t p_reversalPot_Ek = -0.07045
    cdef np.float64_t p_reversalPotCorrection = 0.04
    cdef np.float64_t p_totalCapacitance_C = 6e-12
    cdef np.float64_t p_restingConductance_G0 = 1.974e-09
    cdef np.float64_t p_kConductance_Gk = 1.8e-08
    cdef np.float64_t p_maxMConductance_Gmax = 8e-09
    cdef np.float64_t p_ciliaTimeConst_tc = 0.00213
    cdef np.float64_t p_referencePot = 0.0
    cdef np.float64_t p_sensitivity_s0 = 85e-09
    cdef np.float64_t p_sensitivity_s1 = 5e-09
    cdef np.float64_t p_offset_u0 = 7e-09
    cdef np.float64_t p_offset_u1 = 7e-09

    cdef np.float64_t p_C_ST = 4.0 # IHC bundle compliance in m/N
    cdef np.float64_t p_F0_ST = 1.e3 # IHC bundle stiffness-fluid friction corner frequency to calculate R_ST (Ns/m)

    cdef np.float64_t Z_ST, g_ST, uIHC_old
    cdef np.float64_t dtOverC, gkEpk
    cdef np.float64_t restingPotential_V0

    cdef np.float64_t L_ST, R_ST, f0_ST
    cdef np.float64_t dt

    cdef np.float64_t leakageConductance_Ga,
    cdef np.float64_t conductance_G
    cdef np.float64_t potential_V
    cdef np.float64_t ciliaAct
    cdef np.float64_t u0, u1, s0, s1
    cdef np.float64_t b1_ST, b2_ST, b3_ST
    cdef np.float64_t xST

    cdef Py_ssize_t i

    cdef int sections = 100

    restingPotential_V0 = ((p_restingConductance_G0 *
                            p_endocochlearPot_Et + p_kConductance_Gk *
                            (p_reversalPot_Ek + p_endocochlearPot_Et *
                             p_reversalPotCorrection)) /
                           (p_restingConductance_G0 + p_kConductance_Gk))


    dt = 1 / fs
    dtOverC = dt / p_totalCapacitance_C
    gkEpk = p_kConductance_Gk * (p_reversalPot_Ek +
                                 p_endocochlearPot_Et * p_reversalPotCorrection)

    Z_ST = 0

    L_ST = 1 / p_C_ST            # v-U / F-I analogy !!!
    f0_ST = 2000 * (( (2000/200)**(1/sections) )**(-sec)) # grade Fluid filter from 2000 to 200 Hz

    R_ST = L_ST * 2 * np.pi * f0_ST     # Fluid friction on bundles Ns/m

    g_ST = R_ST / (R_ST + 2*fs*L_ST)  # gamma

    s0 = p_sensitivity_s0
    u0 = p_offset_u0
    s1 = p_sensitivity_s1
    u1 = p_offset_u1


    ciliaAct = 1.0 / (1.0 + exp(u0 / s0) * ( 1 + exp(u1 / s1)))
    leakageConductance_Ga = (p_restingConductance_G0 -
                             p_maxMConductance_Gmax * ciliaAct)

    uIHC_old = restingPotential_V0

    for i in range(len(xbm)):
        b3_ST = -(xbm[i] + Z_ST)
        b1_ST = xbm[i] + 2*g_ST*b3_ST
        b2_ST = -(2*xbm[i] + 2*g_ST*b3_ST + Z_ST)
        xST = -(b2_ST+Z_ST)/2

        Z_ST = -b2_ST


        xST *= ciliaCouplingGain[sec]


        potential_V = uIHC_old

        ciliaAct = 1 / (1 + exp((u0 - xST) / s0) *
                        (1 + exp((u1 - xST) / s1)))

        conductance_G = (p_maxMConductance_Gmax * ciliaAct +
                         leakageConductance_Ga)


        uIHC[i] = (potential_V - dtOverC *
                   (conductance_G * (potential_V - p_endocochlearPot_Et) +
                    p_kConductance_Gk * potential_V - gkEpk))


        uIHC_old = uIHC[i]

    return uIHC
