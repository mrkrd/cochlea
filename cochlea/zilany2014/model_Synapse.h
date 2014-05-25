double Synapse(double *ihcout, double tdres, double cf, int totalstim, int nrep, double spont, double noiseType, double implnt, double sampFreq, double *synouttmp);
int SpikeGenerator(double *synouttmp, double tdres, int totalstim, int nrep, double *sptime);
