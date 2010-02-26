/* Modified version of HCRPSham3StVIn (IHCRP DSAM module) */
#include <stdio.h>
#include "ihcrp.h"


double p_endocochlearPot_Et = 0.1;
double p_reversalPot_Ek = -0.07045;
double p_reversalPotCorrection = 0.04;
double p_totalCapacitance_C = 6e-12;
double p_restingConductance_G0 = 1.974e-09;
double p_kConductance_Gk = 1.8e-08;
double p_maxMConductance_Gmax = 8e-09;
double p_ciliaTimeConst_tc = 0.00213;
/* double p_ciliaCouplingGain_C = 16; /\* [dB] *\/ */
double p_referencePot = 0.0;
double p_sensitivity_s0 = 85e-09;
double p_sensitivity_s1 = 5e-09;
double p_offset_u0 = 7e-09;
double p_offset_u1 = 7e-09;


double p_C_ST = 4.0;            /*  IHC bundle compliance in m/N */
double p_F0_ST = 1.e3;         /*  IHC bundle stiffness-fluid friction corner frequency to calculate R_ST (Ns/m) */

/* x_ST[SECTIONS], */
double Z_ST[SECTIONS], g_ST[SECTIONS];
double *uIHC_old;
/* double *xST_old */
double dtOverC, gkEpk;
double  restingPotential_V0;

/* FILE *debug_fp; */


/* Initialize parameters */
double ihcrp_init(double f_s)
{
     int sec;
     double L_ST,R_ST,f0_ST;
     double dt;

     /* debug_fp = fopen("debug.dat", "w"); */

     restingPotential_V0 = (p_restingConductance_G0 *
			    p_endocochlearPot_Et + p_kConductance_Gk *
			    (p_reversalPot_Ek + p_endocochlearPot_Et *
			     p_reversalPotCorrection)) /
	  (p_restingConductance_G0 + p_kConductance_Gk);

     /* uIHC_old = mxMalloc(SECTIONS * sizeof(double)); */
     /* for (sec = 0; sec < SECTIONS; sec++) { */
     /* 	  /\*** Put reset (to zero ?) code here ***\/ */
     /* 	  uIHC_old[sec]  = restingPotential_V0; */
     /* } */

     dt = 1./f_s;
     dtOverC = dt / p_totalCapacitance_C;
     gkEpk = p_kConductance_Gk * (p_reversalPot_Ek +
				  p_endocochlearPot_Et * p_reversalPotCorrection);

     for(sec=0;sec<SECTIONS;sec++)
     {
	  /* - - - - - - - - - - -  Inner hair cells:  Bundle stimulation  - - - - - - */
	  Z_ST[sec]=0;

	  /*  C_ST=4; ****** */                             /* IHC bundle compliance in m/N */
	  L_ST=1/p_C_ST;                              /* v-U / F-I analogy !!! */
	  f0_ST=2000. * pow(pow(2000./200.,1./SECTIONS),-sec);   /* grade Fluid filter from 2000 to 200 Hz */

	  R_ST=L_ST*2*M_PI*f0_ST;                     /* Fluid friction on bundles Ns/m */

	  g_ST[sec]=R_ST/(R_ST+2*f_s*L_ST);         /* gamma */

     }

     return restingPotential_V0;
}



void
ihcrp(double *uIHC, double *xBM, double *ciliaCouplingGain)
{
     int sec;
     double  leakageConductance_Ga, conductance_G, potential_V;
     /* double  ciliaDisplacement_u; */
     double  ciliaAct;
     double  u0, u1, s0, s1;
     double b1_ST, b2_ST, b3_ST;
     int first_sec = 0;
     int last_sec = SECTIONS;
     double xST[SECTIONS];

     s0 = p_sensitivity_s0;
     u0 = p_offset_u0;
     s1 = p_sensitivity_s1;
     u1 = p_offset_u1;

     ciliaAct = 1.0 / (1.0 + exp(u0 / s0) * ( 1 + exp(u1 / s1)));
     leakageConductance_Ga = p_restingConductance_G0 -
	  p_maxMConductance_Gmax * ciliaAct;

     /* cGain = pow(10.0, p_ciliaCouplingGain_C / 20.0); */
/*    dtOverTc = dt / p_ciliaTimeConst_tc; */



     for(sec=last_sec-1;sec>=first_sec;sec--) {

	  /* - - - - - - - - - - -  Inner hair cells:  Bundle stimulation in fluid UPDATE - - */
	  b3_ST=-(xBM[sec]+Z_ST[sec]);
	  b1_ST=xBM[sec]+2*g_ST[sec]*b3_ST;
	  b2_ST=-(2*xBM[sec]+2*g_ST[sec]*b3_ST+Z_ST[sec]);
	  xST[sec]=-(b2_ST+Z_ST[sec])/2;

	  Z_ST[sec]=-b2_ST;                /* Z storage of bundle highpass filter */


	  /* ciliaDisplacement_u = cGain*xST[sec]; */
	  /*ciliaDisplacement_u = ciliaCouplingGain[sec]*xST[sec]; */
	  xST[sec] *= ciliaCouplingGain[sec];



	  if (uIHC_old == NULL) { /* first run */
	       potential_V = restingPotential_V0;
	  } else {
	       potential_V = uIHC_old[sec];
	  }


/*         ciliaDisplacement_u = xST_old[sec];        only needed with Shamma xST*/

/* Shamma, input is BM velocity! */
/*            ciliaDisplacement_u += dt * cGain * (*inPtr ) -
              ciliaDisplacement_u * dtOverTc; */

	  /* ciliaAct = 1.0 / (1.0 + exp((u0 - ciliaDisplacement_u) / s0) *
	     (1.0 + exp((u1 - ciliaDisplacement_u) / s1))); */
	  ciliaAct = 1.0 / (1.0 + exp((u0 - xST[sec]) / s0) *
			    (1.0 + exp((u1 - xST[sec]) / s1)));

	  conductance_G = p_maxMConductance_Gmax * ciliaAct +
	       leakageConductance_Ga;



	  uIHC[sec] = potential_V - dtOverC *
	       (conductance_G * (potential_V - p_endocochlearPot_Et) +
		p_kConductance_Gk * potential_V - gkEpk);


	  /* if (sec == 10) uIHC[sec] = 100; */
     }
     uIHC_old = uIHC; /* set pointer to just calculated data */

     /* fwrite(uIHC, sizeof(double), SECTIONS, debug_fp); */

     /* xST_old = xST; only needed with Shamma's xST */
}
