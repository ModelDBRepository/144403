/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_LTVNEURON_CC
#define CN_LTVNEURON_CC

#include "CN_neuron.cc"

LTVneuron::LTVneuron(int inlabel, double *the_p= LTV_p):
  neuron(inlabel, LTV_IVARNO, LTVNEURON, the_p, LTV_PNO)
{
}

LTVneuron::LTVneuron(int inlabel, vector<int> inpos, double *the_p= LTV_p):
  neuron(inlabel, LTV_IVARNO, LTVNEURON, inpos, the_p, LTV_PNO)
{
}

inline double LTVneuron::E(double *x)
{
  assert(enabled);
  return x[idx];
}

void LTVneuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }
  
  // differential eqn for E, the membrane potential
  dx[idx]= x[idx]*(1.0 - p[0]*x[idx]) - Isyn;
}

void LTVneuron::noise(double *x, double *dx)
{
  // differential eqn for E, the membrane potential
  dx[idx]= p[1]*RG.n();
}

#endif
