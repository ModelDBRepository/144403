/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_VALNEURON_CC
#define CN_VALNEURON_CC

#include "CN_neuron.cc"

Colpitts::Colpitts(int inlabel, double *the_p= Colp_p):
  neuron(inlabel, Colp_IVARNO, VALNEURON, the_p, Colp_PNO)
{
}

Colpitts::Colpitts(int inlabel, vector<int> inpos, double *the_p= Colp_p):
  neuron(inlabel, Colp_IVARNO, VALNEURON, inpos, the_p, Colp_PNO)
{
}

inline double Colpitts::E(double *x)
{
  assert(enabled);
  return x[idx];
}

void Colpitts::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }

  dx[idx] = p[0]*x[idx+1] + Isyn;
  dx[idx+1] = -p[1]*(x[idx]+x[idx+2]) -p[2]*x[idx+1];
  dx[idx+2] = p[3]*(x[idx+1] + 1.0 - exp(-x[idx]));
}

#endif
