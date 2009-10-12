#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ffGn.h"

int main(void)
{
     int times = 10000;

     /* ffGn args */
     int N = 1000;
     double Hinput = 0.9;
     double mu = 1.0;

     double *randNums;
     int i;


     /* pyResample args */
     int old_len = 300;
     double x[old_len];
     int p = 6;
     int q = 2;
     double *y;


     /* pyRand args */
     int len = 5000;
     double *z;


     printf("*********** ffGn ************\n");
     for (i=0; i<times; i++) {
	  randNums = ffGn(N, Hinput, mu);
	  free(randNums);
     }

     randNums = ffGn(N, Hinput, mu);
     for (i=0; i<N; i++) {
     	  printf("%f\n", randNums[i]);
     }
     free(randNums);


     printf("\n");

     randNums = ffGn(N, Hinput, mu);

     for (i=0; i<N; i++) {
     	  printf("%f\n", randNums[i]);
     }
     free(randNums);

     /* printf("Testing for memory leaks...\n"); */

     /* for (i=0; i<1000000; i++) { */
     /* 	  randNums = ffGn(N*10, Hinput, mu); */
     /* 	  free(randNums); */
     /* } */


     printf("*********** pyResample ************\n");
     for (i=0; i<old_len; i++) {
	  x[i] = i;
     }
     for (i=0; i<times; i++) {
	  y = pyResample(x, old_len, p, q);
	  free(y);
     }

     y = pyResample(x, old_len, p, q);
     for (i=0; i<ceil(old_len*p/q); i++) {
	  printf("%f\n", y[i]);
     }

     printf("\n");

     for (i=0; i<old_len; i++) {
	  x[i] = i;
     }
     y = pyResample(x, old_len, p, q);
     for (i=0; i<ceil(old_len*p/q); i++) {
	  printf("%f\n", y[i]);
     }


     printf("*********** pyRand ************\n");
     for (i=0; i<times; i++) {
	  z = pyRand(len);
	  free(z);
     }

     z = pyRand(len);
     for (i=0; i<len; i++) {
	  printf("%f\n", z[i]);
     }

     return 0;
}
