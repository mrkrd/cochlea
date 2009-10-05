#include <stdio.h>
#include "ffGn.h"

int main(void)
{

     /* ffGn args */
     int N = 10;
     double Hinput = 0.9;
     double mu = 1.0;

     double *randNums;
     int i;


     /* pyResample args */
     int old_len = 3;
     double x[old_len];
     int new_len = 6;
     double *y;


     printf("*********** ffGn ************\n");
     randNums = ffGn(N, Hinput, mu);

     for (i=0; i<N; i++) {
     	  printf("%f\n", randNums[i]);
     }


     randNums = ffGn(N, Hinput, mu);

     for (i=0; i<N; i++) {
     	  printf("%f\n", randNums[i]);
     }


     printf("*********** pyResample ************\n");
     for (i=0; i<old_len; i++) {
	  x[i] = i;
     }
     y = pyResample(x, old_len, new_len);
     for (i=0; i<new_len; i++) {
	  printf("%f\n", y[i]);
     }

     return 0;
}
