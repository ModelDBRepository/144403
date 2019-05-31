/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#include "CN_ECdemiGapsynapse.h"
#include "CN_synapse.cc"

ECdemiGapsynapse:: ECdemiGapsynapse(neuron* source, neuron *target, double ingsyn, double intau, double inVmid, double inVslope):
  synapse(source, target, ECdemiGapIVARNO, ECdemiGapPNO, ECDEMIGAPSYNAPSE)
{
  p[0]= ingsyn;               // injected current
  p[1]= intau;
  p[2]= inVmid;
  p[3]= inVslope;
}


ECdemiGapsynapse::~ECdemiGapsynapse()
{
}

double ECdemiGapsynapse::gsyn()
{
  return p[0];
}

void ECdemiGapsynapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double ECdemiGapsynapse::Isyn(double *x)
{
  return p[0]*x[idx]*(source->E(x) - target->E(x));
}

void ECdemiGapsynapse::derivative(double *x, double *dx)
{
  static double Sinf;
  Sinf= 1.0/(1.0+exp((target->E(x)-p[2])/p[3]));
  dx[idx]= -(x[idx]-Sinf)/p[1];
}



