/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_VDPOLNEURON_CC
#define CN_VDPOLNEURON_CC

#include "CN_neuron.cc"

VdPolneuron::VdPolneuron(int inlabel, double *the_p= VDPOL_p):
  neuron(inlabel, VDPOL_IVARNO, VDPOLNEURON, the_p, VDPOL_PNO)
{
}

VdPolneuron::VdPolneuron(int inlabel, vector<int> inpos, double *the_p= VDPOL_p):
  neuron(inlabel, VDPOL_IVARNO, VDPOLNEURON, inpos, the_p, VDPOL_PNO)
{
}

inline double VdPolneuron::E(double *x)
{
  assert(enabled);
  return x[idx];
}

void VdPolneuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }
  
  // differential eqn for E, the membrane potential
  dx[idx]= x[idx+1];
  dx[idx+1]= (p[0]-x[idx]*x[idx])*x[idx+1] - p[1]*x[idx] + Isyn;
}

void VdPolneuron::noise(double *x, double *dx)
{
  // differential eqn for E, the membrane potential
  dx[idx]= p[2]*RG.n();
}

#endif
