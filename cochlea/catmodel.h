
void IHCAN(double *px,
           double cf,
           int nrep,
           double tdres,
           int totalstim,
           double cohc,
           double cihc,
           double *ihcout);

double Synapse(double *ihcout,
	       double tdres,
	       double cf,
	       int totalstim,
	       int nrep,
	       double spont,
	       double implnt,
	       double sampFreq,
	       double *synouttmp,
	       int with_ffGn);

int SpikeGenerator(double *synouttmp,
		   double tdres,
		   int totalstim,
		   int nrep,
		   double *sptime);
