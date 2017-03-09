
/********  header file, arrays and functions common to all following subroutines and
  the same for all problems.  NOTE this file preceeded by compiler directives defining
  parameters and the track*[] arrays
********/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>

/* global arrays for history and geometry */

static double history[SIZE][NSTEP][NCELLTOT];
static int geometry[NCELLTOT][NNEIGHBOR];

double MAX(double a,double b){
	 if (a>b) return a;
 	 return b;
}
 
double FRAND()  {
	return (double) rand()/((double)RAND_MAX + 1);
}
 
double MIN( double a, double b ) {
	if (a<b) return a;
	return b;
}



double POW(double x,double n){
  return exp(n*log(x));
}





double HillR(double x,double thresh,double n)
{
	double r=exp(n*log(x/thresh));
	return 1.0/(1+r);
}
 
 
double HillA(double x,double thresh,double n)
{
	double r=exp(n*log(x/thresh));
	return r/(1+r);
 }

