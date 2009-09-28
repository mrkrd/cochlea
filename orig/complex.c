/** 
complex.cpp includes all of the COMPLEX math functions needed for model programs
*/

#include <stdlib.h>
#include <math.h>
#include "complex.hpp"

//divide
COMPLEX compdiv(COMPLEX ne,COMPLEX de)
{
  double d;
  COMPLEX z;
  d=de.x*de.x+de.y*de.y;
  z.x=(ne.x*de.x+ne.y*de.y)/d;
  z.y=(ne.y*de.x-ne.x*de.y)/d;
  return(z);
}
// this returns a complex number equal to exp(i*theta)
COMPLEX compexp(double theta)
{
  COMPLEX dummy;
  dummy.x = cos(theta);
  dummy.y = sin(theta);
  return dummy;
}
// Multiply a complex number by a scalar
COMPLEX compmult(double scalar, COMPLEX compnum)
{
 COMPLEX answer;
 answer.x = scalar * compnum.x;
 answer.y = scalar * compnum.y;
 return answer;
}
// Find the product of 2 complex numbers
COMPLEX compprod(COMPLEX compnum1, COMPLEX compnum2)
{
 COMPLEX answer;
 answer.x = (compnum1.x * compnum2.x) - (compnum1.y * compnum2.y);
 answer.y = (compnum1.x * compnum2.y) + (compnum1.y * compnum2.x);
 return answer;
}
// add 2 complex numbers
COMPLEX comp2sum(COMPLEX summand1, COMPLEX summand2)
{
 COMPLEX answer;
 answer.x = summand1.x + summand2.x;
 answer.y = summand1.y + summand2.y;
 return answer;
}
// add three complex numbers
COMPLEX comp3sum(COMPLEX summand1, COMPLEX summand2, COMPLEX summand3)
{
 COMPLEX answer;
 answer.x = summand1.x + summand2.x + summand3.x;
 answer.y = summand1.y + summand2.y + summand3.y;
 return answer;
}

// subtraction: complexA - complexB
COMPLEX compsubtract(COMPLEX complexA, COMPLEX complexB)
{
 COMPLEX answer;
 answer.x = complexA.x - complexB.x;
 answer.y = complexA.y - complexB.y;
 return answer;
}
//Get the real part of the complex
double REAL(COMPLEX compnum)
{ return(compnum.x); };

//Get the imaginary part of the complex
double IMAG(COMPLEX compnum)
{ return(compnum.y); };

//Get the conjugate of the complex signal
COMPLEX compconj(COMPLEX complexA)
{
  COMPLEX answer;
  answer.x = complexA.x;
  answer.y = -complexA.y;
  return (answer);
};