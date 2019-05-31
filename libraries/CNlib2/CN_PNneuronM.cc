/*--------------------------------------------------------------------------
   Author: Thomas Nowotny

   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402

   email to:  tnowotny@ucsd.edu

   initial version: 2005-08-17

--------------------------------------------------------------------------*/

#ifndef CN_PNNEURON_CC
#define CN_PNNEURON_CC

#include "CN_neuron.cc"

PNneuron::PNneuron(int inlabel, double *the_p= PN_p):
  neuron(inlabel, PN_IVARNO, PNNEURON, the_p, PN_PNO)
{
}

PNneuron::PNneuron(int inlabel, vector<int> inpos, double *the_p= PN_p):
  neuron(inlabel, PN_IVARNO, PNNEURON, inpos, the_p, PN_PNO)
{
}

inline double PNneuron::E(double *x)
{
  assert(enabled);
  return x[idx];
}


void PNneuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }

  // differential eqn for E, the membrane potential
  dx[idx]= -(pw3(x[idx+1])*x[idx+2]*p[0]*(x[idx]-p[1]) +
	     pw4(x[idx+3])*p[2]*(x[idx]-p[3])+
	     p[4]*(x[idx]-p[5])+p[6]*(x[idx]-p[7])+
	     p[10]*x[idx+4]*(x[idx]-p[3])-Isyn-p[11])/p[9];

  // diferential eqn for m, the probability for one Na channel activation
  // particle
  _a= 0.32*(13.0-x[idx]-p[8]) / (exp((13.0-x[idx]-p[8])/4.0)-1.0);
  _b= 0.28*(x[idx]+p[8]-40.0)/(exp((x[idx]+p[8]-40.0)/5.0)-1.0);
  dx[idx+1]= _a*(1.0-x[idx+1])-_b*x[idx+1];

  // differential eqn for h, the probability for the Na channel blocking
  // particle to be absent
  _a= 0.128*exp((17.0-x[idx]-p[8])/18.0);
  _b= 4.0 / (exp((40-x[idx]-p[8])/5.0)+1.0);
  dx[idx+2]= _a*(1.0-x[idx+2])-_b*x[idx+2];

  // differential eqn for n, the probability for one K channel activation
  // particle
  _a= .032*(15.0-x[idx]-p[8]) / (exp((15.0-x[idx]-p[8])/5.0)-1.0);
  _b= 0.5*exp((10.0-x[idx]-p[8])/40.0);
  dx[idx+3]= _a*(1.0-x[idx+3])-_b*x[idx+3];

  // M current activation
  dx[idx+4]= tauz*(1.0/(1.0+exp(-(x[idx]+20.0)/5.0)) - x[idx+4]);
}


#endif
