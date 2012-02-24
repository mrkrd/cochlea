#define SECTIONS 100

void LCR4_init_c(double f_s,
		 double *freq_map,
		 double *Qmin,
		 double *SAT1,
		 double *SAT4);


void LCR4_c(double *xBM,
	    double *Qmax,
	    double *Qmin,
	    int first_sec,
	    int last_sec);
