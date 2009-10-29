#ifndef _COMPLEX_HPP
#define _COMPLEX_HPP

/* COMPLEX.HPP header file		
 * use for complex arithmetic in C 
 (part of them are from "C Tools for Scientists and Engineers" by L. Baker)
*/

///Structure of the complex
struct __COMPLEX{ double x,y; };
/// structure COMPLEX same as __COMPLEX
typedef struct __COMPLEX COMPLEX;

/* for below, X, Y are complex structures, and one is returned */

///real part of the complex multiplication
#define CMULTR(X,Y) ((X).x*(Y).x-(X).y*(Y).y)
///image part of the complex multiplication
#define CMULTI(X,Y) ((X).y*(Y).x+(X).x*(Y).y)
/// used in the Division : real part of the division
#define CDRN(X,Y) ((X).x*(Y).x+(Y).y*(X).y)
/// used in the Division : image part of the division
#define CDIN(X,Y) ((X).y*(Y).x-(X).x*(Y).y)
/// used in the Division : denumerator of the division
#define CNORM(X) ((X).x*(X).x+(X).y*(X).y)
///real part of the complex
#define CREAL(X) (double((X).x))
///conjunction value
#define CONJG(z,X) {(z).x=(X).x;(z).y= -(X).y;}
///conjunction value
#define CONJ(X) {(X).y= -(X).y;}
///muliply : z could not be same variable as X or Y, same rule for other Macro
#define CMULT(z,X,Y) {(z).x=CMULTR((X),(Y));(z).y=CMULTI((X),(Y));}
///division
#define CDIV(z,X,Y){double d=CNORM(Y); (z).x=CDRN(X,Y)/d; (z).y=CDIN(X,Y)/d;}
///addition
#define CADD(z,X,Y) {(z).x=(X).x+(Y).x;(z).y=(X).y+(Y).y;}
///subtraction
#define CSUB(z,X,Y) {(z).x=(X).x-(Y).x;(z).y=(X).y-(Y).y;}
///assign
#define CLET(to,from) {(to).x=(from).x;(to).y=(from).y;}
///abstract value(magnitude)
#define cabs(X) sqrt((X).y*(X).y+(X).x*(X).x)
///real to complex
#define CMPLX(X,real,imag) {(X).x=(real);(X).y=(imag);}
///multiply with real
#define CTREAL(z,X,real) {(z).x=(X).x*(real);(z).y=(X).y*(real);}


/* implementation using function : for compatibility */
///divide
COMPLEX compdiv(COMPLEX ne,COMPLEX de);
///this returns a complex number equal to exp(i*theta)
COMPLEX compexp(double theta);
/// Multiply a complex number by a scalar
COMPLEX compmult(double scalar,COMPLEX compnum);
/// Find the product of 2 complex numbers
COMPLEX compprod(COMPLEX compnum1, COMPLEX compnum2);
/// add 2 complex numbers
COMPLEX comp2sum(COMPLEX summand1, COMPLEX summand2);
/// add three complex numbers
COMPLEX comp3sum(COMPLEX summand1, COMPLEX summand2, COMPLEX summand3);
/// subtraction: complexA - complexB
COMPLEX compsubtract(COMPLEX complexA, COMPLEX complexB);
///Get the real part of the complex
double  REAL(COMPLEX compnum);//{double(compnum.x);};
//Get the imaginary part of the complex
double IMAG(COMPLEX compnum);
/// Get the conjugate of the complex 
COMPLEX compconj(COMPLEX complexA);


#endif
