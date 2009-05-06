/*  ----------------------------------------------------------------------------- */
/*  ------- ¨2000 Werner Hemmert, University of Z’rich and IBM ZRL -------------- */
/*  ----------------------------------------------------------------------------- */
/*  ------- modified ¨2001 Werner Hemmert, Infineon Corporate Research ---------- */
/*  ----------------------------------------------------------------------------- */
/*  ----------------------------------------------------------------------------- */
/*  ------- modified 2002-3 Marcus Holmberg, Infineon Corporate Research ---------- */
/*  ----------------------------------------------------------------------------- */
/*  Wave Digital Filter of BM (based on bm_Strube2.m) */
/*    Broken up into subroutines to be linked and executed from model.c */
/*  7.6.2000 Converted to "C" */
/*  18.6.2000 increase Q with frequency */
/*  31.8.2000 Broken up into subroutines to be linked and executed from model.c */
/*  23.05.2001 C++  output is xBM */
/*  17.10.2002 C convention for indexing (0..SECTIONS-1) */
/*  11.12.2002 freq_map and Q as input parameter */
/*  12.12.2002 algorithm improved, less memory usage */
/*  31.03.2003 remove tanh correction (CORRECTION SHIFTS FREQU RANGE TOO FAR AWAY!) and increase L_S (LARGER L_S IMPROVES HF-SLOPE) */
/*  08.05.2003 new parameter values, derived from Viergever (Ph.D. thesis, 1980) */
/* 15.05.2003 Changed to work with nonlinear model */
/* 15.05.2003 Bug in calculation of port resistances corrected */
/* 26.05.2003 Parameter 'mode' introduced */
/*            mode = 0: standard */
/*            mode = 1: output is linear BM */
/* 28.05.2003 Program divided into "quick" and "debug" (depending on 'debug_on') */
/* 26.8.2003  Uses whole basal end of the cochlea */
/* 30.9.2003  Uses frequencies up to 15kHz */
/*  3.11.2003 New program structure with parameter file */
/*  2.5.2004  Uses an ideal current source as input */

/* Global variables for calculation */
double g11[SECTIONS], g12[SECTIONS], g2[SECTIONS], g3[SECTIONS], g41[SECTIONS], g42[SECTIONS], gh;
double Z13[SECTIONS], Z32[SECTIONS], Z42[SECTIONS], Z43[SECTIONS], Zhel;
/* This is a fix to be able to calculate input impedance of the model */
double R_input;

/*--------------------------------------- INIT ------------------------------*/
double bm_init(double f_s, double *Ls, double *Rs, double *Ct, double *Rbm, double *Cbm, double *Lbm, double Rh, double Lh)
{
     int sec;
     double R14, R44, G33, G23;

     for(sec=SECTIONS-1;sec>=0;sec--)
     {

	  if (sec == (SECTIONS-1)) {
	       /* HELICOTREMA */
	       R14 = Rh + (2.*f_s*Lh);
	       gh = Rh/R14;
	  }

	  /*------------------------------------adaptor 4 (series)---------------------------------------*/
	  R44 = Rbm[sec] + 2*f_s*Lbm[sec] + 1./(2*f_s*Cbm[sec]);
	  g41[sec] = Rbm[sec]/R44;
	  g42[sec] = 2*f_s*Lbm[sec]/R44;

	  /*----------------------------------adaptor 3 (parallel)--------------------------------------*/
	  G33 = 1/R44 + 2*f_s*Ct[sec];
	  g3[sec] = 1./(G33*R44);

	  /*----------------------------------adaptor 2 (parallel)--------------------------------------*/
	  G23 = 1/R14 + G33;
	  g2[sec] = 1./(R14*G23);

	  /*---------------------------------adaptor 1 (series)--------------------------------------------*/
	  R14 = 1./G23 + Rs[sec] + 2*f_s*Ls[sec];
	  g11[sec] = 1./(G23*R14);
	  g12[sec] = Rs[sec]/R14;

	  /*---------Initialize Storage elements-----------*/
	  Z13[sec]=0;
	  Z32[sec]=0;
	  Z42[sec]=0;
	  Z43[sec]=0;
     }
     Zhel = 0.;
     R_input = R14;
     /*------- mexPrintf("R_input = %f\n",R_input);---------*/
     return 1;
}

/*---------------------------------------------------------------------------*/
/*------------------------------- E N D - INIT ------------------------------*/
/*---------------------------------------------------------------------------*/

#if debug_on
void bm_wave(double input, double *x_BM, double *debug, double *ampl_corr, double *Abm, double *Cbm,int mode)
#endif
#if !debug_on
     void bm_wave(double input, double *x_BM, double *ampl_corr, double *Abm, double *Cbm)
#endif
{
     int sec;
     double  b14[SECTIONS], b44[SECTIONS],b30[SECTIONS], b33[SECTIONS], b20[SECTIONS], b23[SECTIONS];
     double a10, b11, b12, b13,
	  a40, b41, b42, b43,
	  b21, b22,
	  b31, b32,
	  ah0,bh1,bh2,
	  a21[SECTIONS], a14;

/*-------------------------------------- backward wave ----------------------*/


     for(sec=SECTIONS-1;sec>=0;sec--) {
	  if(sec==(SECTIONS-1)) {
	       /* a21[sec]=Zhel; */
	       a21[sec]=Zhel;
	  } else {
	       /* a21[sec]=b14[sec+1]; */
	       a21[sec]=-b14[sec+1];
	  }

	  b44[sec] = -(-Z42[sec] + Z43[sec]);

	  b30[sec] = -g3[sec]*(Z32[sec] - b44[sec]);
	  b33[sec]  = Z32[sec] + b30[sec];

	  b20[sec] = -g2[sec]*(b33[sec]-a21[sec]);
	  b23[sec] = b33[sec] + b20[sec];

	  b14[sec] = -(b23[sec] - Z13[sec]);

     }
     /*-------------------------------------- forward wave ----------------------*/

     for(sec=0;sec<SECTIONS;sec++) {

	  if (sec == 0) {
	       a14 = 2.*R_input*input + b14[0];
#if debug_on
	       if (mode == 4) {
		    /* Debug[0,:] contains  voltage over input */
		    debug[0] = (a14 + b14[0])/2.;
		    /* Debug[1,:] contains input current */
		    debug[1] = (a14-b14[0])/(2.*R_input);
	       }
#endif
	       /* mexPrintf("Input signal = %f\n",input);
		  mexPrintf("Input wave = %f\n",a14); */
	  } else {
	       /* a14 = b21; */
	       a14 = -b21;

	  }

	  a10 = a14 - b14[sec];
	  b11 = b23[sec] - g11[sec]*a10;
	  b12 = -g12[sec]*a10;
	  b13 = -(b11 + b12 + a14);

	  b22 = b11 + b20[sec];
	  b21 = b22 + b33[sec] - a21[sec];

	  b22      = b11 + b20[sec];
	  b21      = b22 + b33[sec] - a21[sec];

	  b32      = b22 + b30[sec];
	  b31      = b32 + Z32[sec] - b44[sec];

	  a40      = b31 - b44[sec];
	  b41      = -g41[sec]*a40;
	  b42      = -Z42[sec] - g42[sec]*a40;
	  b43      = -(b41 + b42 + b31);

	  /* v_BM[sec-1] =(b22[sec]-b34[sec])/2./Z[sec];               BM Velocity */
	  /* x_BM[sec-1]+=(b22[sec]-b34[sec])/2./Z[sec]/f_s;           BM dislpacement */
	  /* x_BM[sec-1]=(b33[sec]+Z3[sec])/4.*C[sec];                 BM dislpacement x = F * C */

	  x_BM[sec] = (b43+Z43[sec])*Cbm[sec]/2./Abm[sec]*ampl_corr[sec];


#if debug_on
	  if (mode == 1 || mode == 9) {
	       debug[sec] = (b43+Z43[sec])*Cbm[sec]/2./Abm[sec]; /* return linear BM without ampl_corr, "abort" calculation (xBM = 0) */
	  }
#endif

	  Z13[sec] = b13;
	  Z32[sec] = b32;
	  Z42[sec] = b42;
	  Z43[sec] = b43;

     }
     /* - - - - - - - Helicotrema - - - - - - - - - */
     ah0 = b21 - Zhel;
     bh1 = -gh * ah0;
     bh2 = -(bh1 + b21);
     Zhel = bh2;
}

/*---------------------------------------------------------------------------*/
/*--------------------------------- BM - E N D ------------------------------*/
/*---------------------------------------------------------------------------*/
