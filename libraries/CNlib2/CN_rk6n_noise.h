/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------

  Implementation of a 6-5 Runge Kutta method with adaptive time step
  mostly taken from the book "The numerical analysis of ordinary differential
  equations - Runge-Kutta and general linear methods" by J.C. Butcher, Wiley,
  Chichester, 1987 and a free adaption to a 6 order Runge Kutta method
  of an ODE system with additive white noise

--------------------------------------------------------------------------*/  

using namespace std;

#ifndef CN_RK6N_NOISE_H
#define CN_RK6N_NOISE_H

#include "CN_NeuronModel.h"
#include <cmath>

class rk6n_noise
{
private:
  double a[9][8];
  double b[9];

  double sqrtdt;
  double *Y[9];
  double *F[9];
  double *NS[9];
  double aF;
  double aNS;
  int i, j, k;
  
protected:
  int N;

public:
  rk6n_noise(int);
  ~rk6n_noise();
  void integrate(double *, double *, NeuronModel *, double);
};

#endif
