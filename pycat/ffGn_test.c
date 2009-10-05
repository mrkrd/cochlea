#include <stdio.h>
#include "ffGn.h"

int main(void)
{

     /* ffGn args */
     int N = 1000;
     double Hinput = 0.9;
     double mu = 1.0;

     double *randNums;
     int i;


     /* pyResample args */
     int old_len = 3;
     double x[old_len];
     int p = 3;
     int q = 2;
     double *y;


     /* pyRand args */
     int len = 5;
     double *z;


     printf("*********** ffGn ************\n");
     randNums = ffGn(N, Hinput, mu);

     for (i=0; i<N; i++) {
     	  printf("%f\n", randNums[i]);
     }


     printf("\n");

     randNums = ffGn(N, Hinput, mu);

     for (i=0; i<N; i++) {
     	  printf("%f\n", randNums[i]);
     }


     printf("*********** pyResample ************\n");
     for (i=0; i<old_len; i++) {
	  x[i] = i;
     }
     y = pyResample(x, old_len, p, q);
     for (i=0; i<ceil(old_len*p/q); i++) {
	  printf("%f\n", y[i]);
     }


     printf("*********** pyRand ************\n");
     z = pyRand(len);
     for (i=0; i<len; i++) {
	  printf("%f\n", z[i]);
     }

     return 0;
}
