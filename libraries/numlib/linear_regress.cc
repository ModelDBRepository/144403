/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Center for Computational Neuroscience and Robotics
              University of Sussex
	     Falmer, Brighton BN1 9QJ, UK 
  
   email to:  T.Nowotny@sussex.ac.uk
  
   initial version: 2007-02-22
  
--------------------------------------------------------------------------*/


#ifndef LINEAR_REGRESS_CC
#define LINEAR_REGRESS_CC

#include <cmath>

double linear_regress_slope(double *x, double ***y, int j, int n)
{
  double xsum= 0.0, ysum= 0.0, xysum= 0.0, xsqsum= 0.0;
  double slope;

  if (n > 1) {
    for (int i= 0; i < n; i++) {
      xsum+= x[i];
      ysum+= y[i][j][0];
      xysum+= x[i]*y[i][j][0];
      xsqsum+= x[i]*x[i];
    }
    slope= (n*xysum-xsum*ysum)/(n*xsqsum - xsum*xsum);
    return slope;
  }
  return 0.0;
}

#endif
