#define SECTIONS 100

double bm_init_c(double f_s,
		 double *Ls,
		 double *Rs,
		 double *Ct,
		 double *Rbm,
		 double *Cbm,
		 double *Lbm,
		 double Rh,
		 double Lh);

void bm_wave_c(double input,
	     double *x_BM,
	     double *ampl_corr,
	     double *Abm,
	     double *Cbm);


