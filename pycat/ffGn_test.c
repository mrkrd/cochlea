#include <stdio.h>
#include "ffGn.h"

int main(void)
{
     int N = 10;
     double Hinput = 0.9;
     double mu = 1.0;

     double *randNums;
     int i;

     randNums = ffGn(N, Hinput, mu);

     for (i=0; i<N; i++) {
	  printf("%f\n", randNums[i]);
     }

     return 0;
}
