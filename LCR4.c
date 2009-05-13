/* 12.11.2002  Modified to work with Qmax as input parameter */
/*  3.12.2002  Modified to work like proposed in patent application */
/* 12.12.2002  Modified to wotk together with model_bm_freq_map.c */
/* 15.05.2003  Modified to work with separate freq_map and element values */
/* 20.05.2003  Bug in resonator 3 fixed */
/* 28.05.2003  Program splitted into "quick" and "debug" */
/* 28.05.2003  Resonator amplification adjusted (ampl = 1 for high amplitudes) */
/* 23.06.2003  Bug in resonator 4 corrected */
/* 26.6.2003   Q stars with Qmin to suppress onset */
/* 12.10.2003  Changed SAT-values to create true HRS-fibers (compression sets in later) */
/* 16.10.2003  Changed back to start compression at approx 1 nm */
/*  3.11.2003 New program structure with parameter file */
/*  4. 2.2004  SAT-values are set separately for each section. */
/* 21. 9.2005  Setting abs to the xBM (i.e. right place) */

#include <math.h>

#include "LCR4.h"

double T4C[SECTIONS], T4L[SECTIONS], T5C[SECTIONS], T5L[SECTIONS], T6C[SECTIONS], T6L[SECTIONS],
     T7C[SECTIONS], T7L[SECTIONS], Rtest[SECTIONS], Rc_res[SECTIONS], Rl_res[SECTIONS],
     Qd[4][SECTIONS],
     Z_CLP1[SECTIONS],Z_CLP2[SECTIONS],Z_CLP3[SECTIONS],Z_CLP4[SECTIONS],g_CLP,
     p_SAT1[SECTIONS],p_SAT2[SECTIONS],p_SAT3[SECTIONS],p_SAT4[SECTIONS];
double *weight_bm;
int weight_bm_n;
int g_firstsec[SECTIONS],g_lastsec[SECTIONS];


void LCR4_init(double f_s,
	       double *freq_map,
	       double *Qmin,
	       double *SAT1,
	       double *SAT4)
{
     int sec, i;
     double R1, R2, tau;

     tau = 1./(2*M_PI*0.8e3);

     R2 = 1./(2*f_s*tau);
     R1 = 1.;

     g_CLP = (R2-R1)/(R2+R1);

     for (sec=0;sec<SECTIONS;sec++) {
	  Rc_res[sec] = Qmin[sec]/(2.*f_s);
	  Rl_res[sec] = Rc_res[sec]*pow(f_s/(M_PI*freq_map[sec]),2);

	  T4C[sec]=0;
	  T4L[sec]=0;
	  T5C[sec]=0;
	  T5L[sec]=0;
	  T6C[sec]=0;
	  T6L[sec]=0;
	  T7C[sec]=0;
	  T7L[sec]=0;

	  for (i=0;i<4;i++) {
	       Qd[i][sec]=Qmin[sec];
	  }

	  p_SAT1[sec] = SAT1[sec];
	  p_SAT4[sec] = SAT4[sec];
	  /* p_SAT2[sec] = pow(10.,log10(p_SAT4[sec]) + log10(p_SAT1[sec]/p_SAT4[sec])/3.);
	     p_SAT3[sec] = pow(10.,log10(p_SAT4[sec]) + log10(p_SAT1[sec]/p_SAT4[sec])*2./3.); */

	  p_SAT2[sec] = pow(10.,log10(p_SAT4[sec]) + log10(p_SAT1[sec]/p_SAT4[sec])*2./3.);
	  p_SAT3[sec] = pow(10.,log10(p_SAT4[sec]) + log10(p_SAT1[sec]/p_SAT4[sec])/3.);

	  Z_CLP1[sec] = 0.;
	  Z_CLP2[sec] = 0.;
	  Z_CLP3[sec] = 0.;
	  Z_CLP4[sec] = 0.;

     }
}




void LCR4(double *xBM,
	  double *Qmax,
	  double *Qmin,
	  int first_sec,
	  int last_sec)
{
     double V2[SECTIONS], V3[SECTIONS], V4[SECTIONS], V5[SECTIONS],
	  X2[SECTIONS], X3[SECTIONS], X4[SECTIONS], X5[SECTIONS],
	  Rr11, Rr13, Gr21, Gr23, gr1, gr2,
	  b11, b12, b13, b20, b21, b22, b23, a0;

     int sec;
     double rect, rect_LP, b2;
     for (sec=first_sec;sec < last_sec; sec++) {


	  Rr11=sqrt(Rl_res[sec]*Rc_res[sec])/Qd[0][sec];      /* sqrt(L/C)/Q */
	  Rr13 = Rr11 + Rl_res[sec];

	  Gr21 = 1./Rr13;
	  Gr23 = Gr21 + 1./Rc_res[sec];

	  gr1=Rr11/Rr13;
	  gr2=Gr21/Gr23;


/*------------------------- Filtering ---------------------------------*/

	  b13=-(xBM[sec]+T4L[sec]);
	  b20=-gr2*(T4C[sec]-b13);
	  b23=b20+T4C[sec];

	  b22=b20+b23;
	  b21=b22+T4C[sec]-b13;
	  a0=b21-b13;
	  b11=xBM[sec]-gr1*a0;
	  b12=-(b11+b21);

	  X2[sec] = (b22+T4C[sec])/2.;
	  V2[sec] = (T4C[sec]-b22)/(2.*Rc_res[sec]);

	  T4L[sec]=-b12;
	  T4C[sec]=b22;

	  rect = 1.- (2. / (1.+exp( -p_SAT1[sec]*fabs(X2[sec])) ) -1.);
	  b2 = rect*(1. + g_CLP) - Z_CLP1[sec]*g_CLP;
	  rect_LP = (b2 + Z_CLP1[sec])/2.;
	  Z_CLP1[sec] = b2;

	  Qd[0][sec] = (Qmax[sec]-Qmin[sec])*rect_LP+Qmin[sec];

/* RESONATOR 2 */
	  /*------------- Calculating port values from old Qd-value ------------------*/
	  Rr11=sqrt(Rl_res[sec]*Rc_res[sec])/Qd[1][sec];      /*sqrt(L/C)/Q */
	  Rr13 = Rr11 + Rl_res[sec];

	  Gr21 = 1./Rr13;
	  Gr23 = Gr21 + 1./Rc_res[sec];

	  gr1=Rr11/Rr13;
	  gr2=Gr21/Gr23;

	  /*------------------------- Filtering ---------------------------------*/
	  b13=-(X2[sec]+T5L[sec]);
	  b20=-gr2*(T5C[sec]-b13);
	  b23=b20+T5C[sec];

	  b22=b20+b23;
	  b21=b22+T5C[sec]-b13;
	  a0=b21-b13;
	  b11=X2[sec]-gr1*a0;
	  b12=-(b11+b21);

	  X3[sec] = (b22+T5C[sec])/2.;
	  V3[sec] = (T5C[sec]-b22)/(2.*Rc_res[sec]);

	  T5L[sec]=-b12;
	  T5C[sec]=b22;

	  rect = 1.- (2. / (1.+exp( -p_SAT2[sec]*fabs(X3[sec])) ) -1.);
	  b2 = rect*(1. + g_CLP) - Z_CLP2[sec]*g_CLP;
	  rect_LP = (b2 + Z_CLP2[sec])/2.;
	  Z_CLP2[sec] = b2;

	  Qd[1][sec]=(Qmax[sec]-Qmin[sec])*rect_LP+Qmin[sec];

/* RESONATOR 3 */
	  /*------------- Calculating port values from old Qd-value ------------------*/
	  Rr11=sqrt(Rl_res[sec]*Rc_res[sec])/Qd[2][sec];      /*sqrt(L/C)/Q */
	  Rr13 = Rr11 + Rl_res[sec];

	  Gr21 = 1./Rr13;
	  Gr23 = Gr21 + 1./Rc_res[sec];

	  gr1=Rr11/Rr13;
	  gr2=Gr21/Gr23;

	  /*------------------------- Filtering ---------------------------------*/

	  b13=-(X3[sec]+T6L[sec]);
	  b20=-gr2*(T6C[sec]-b13);
	  b23=b20+T6C[sec];

	  b22=b20+b23;
	  b21=b22+T6C[sec]-b13;
	  a0=b21-b13;
	  b11=X3[sec]-gr1*a0;
	  b12=-(b11+b21);

	  X4[sec] = (b22+T6C[sec])/2.;
	  V4[sec] = (T6C[sec]-b22)/(2.*Rc_res[sec]);


	  T6L[sec]=-b12;
	  T6C[sec]=b22;

	  rect = 1.- (2. / (1.+exp( -p_SAT3[sec]*fabs(X4[sec])) ) -1.);
	  b2 = rect*(1. + g_CLP) - Z_CLP3[sec]*g_CLP;
	  rect_LP = (b2 + Z_CLP3[sec])/2.;
	  Z_CLP3[sec] = b2;

	  Qd[2][sec]=(Qmax[sec]-Qmin[sec])*rect_LP+Qmin[sec];

/* RESONATOR 4 */
	  /*------------- Calculating port values from old Qd-value ------------------*/
	  Rr11=sqrt(Rl_res[sec]*Rc_res[sec])/Qd[3][sec];      /*sqrt(L/C)/Q */
	  Rr13 = Rr11 + Rl_res[sec];

	  Gr21 = 1./Rr13;
	  Gr23 = Gr21 + 1./Rc_res[sec];

	  gr1=Rr11/Rr13;
	  gr2=Gr21/Gr23;

	  /*------------------------- Filtering ---------------------------------*/

	  b13=-(X4[sec]+T7L[sec]);
	  b20=-gr2*(T7C[sec]-b13);
	  b23=b20+T7C[sec];

	  b22=b20+b23;
	  b21=b22+T7C[sec]-b13;
	  a0=b21-b13;
	  b11=X4[sec]-gr1*a0;
	  b12=-(b11+b21);

	  X5[sec] = (b22+T7C[sec])/2.;
	  V5[sec] = (T7C[sec]-b22)/(2.*Rc_res[sec]);

	  T7L[sec]=-b12;
	  T7C[sec]=b22;

	  rect = 1.- (2. / (1.+exp( -p_SAT4[sec]*fabs(X5[sec])) ) -1.);
	  b2 = rect*(1. + g_CLP) - Z_CLP4[sec]*g_CLP;
	  rect_LP = (b2 + Z_CLP4[sec])/2.;
	  Z_CLP4[sec] = b2;

	  Qd[3][sec]=(Qmax[sec]-Qmin[sec])*rect_LP+Qmin[sec];



/*    *xBM = V5*0.001;*/        /* displacement of last resonator (in meter)*/

	  xBM[sec] = X5[sec];

     } /* loop over sections */

}
