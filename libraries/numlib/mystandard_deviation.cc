/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-03-19
  
--------------------------------------------------------------------------*/

#ifndef STANDARD_DEVIATION_CC
#define STANDARD_DEVIATION_CC

#include <cmath>

double mystandard_deviation(double sum, double sum_sq, int n)
{
  double tmp;
  if (n > 1) {
    tmp= (sum_sq - sum*sum/n)/(n-1);
    if (tmp > 0.0) {
      tmp= sqrt(tmp);
      if (isnan(tmp)) {
	cerr << "sum_sq" << sum_sq << " sum " << sum << " n " << n << endl;
	exit(1);
      }
      return tmp;
    }
  }
  return 0.0;
}

#endif
