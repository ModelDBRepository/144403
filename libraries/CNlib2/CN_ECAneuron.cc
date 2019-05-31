/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_ECNEURON_CC
#define CN_ECNEURON_CC

#include "CN_neuron.cc"

ECneuron::ECneuron(int inlabel, double *the_p= ECN_p):
  neuron(inlabel, ECN_IVARNO, ECNEURON, the_p, ECN_PNO)
{
}

ECneuron::ECneuron(int inlabel, vector<int> inpos, double *the_p= ECN_p):
  neuron(inlabel, ECN_IVARNO, ECNEURON, inpos, the_p, ECN_PNO)
{
}

inline double ECneuron::E(double *x)
{
  assert(enabled);
  return x[idx];
}

void ECneuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }
  
  // differential eqn for E, the membrane potential
  dx[idx]= -((pw3(x[idx+1])*x[idx+2]*p[0]+p[2]*x[idx+4])*(x[idx]-p[1]) +
	     pw4(x[idx+3])*p[3]*(x[idx]-p[4])+
	     p[7]*(x[idx+5]*0.65+x[idx+6]*0.35)*(x[idx]-p[8])+
	     p[5]*(x[idx]-p[6])-Isyn+2.85)/p[9];

  // diferential eqn for m, the probability for one Na channel activation
  // particle
  _a= -0.1*(x[idx]+23)/(exp(-0.1*(x[idx]+23))-1);
  _b= 4*exp(-(x[idx]+48)/18);
  dx[idx+1]= _a*(1.0-x[idx+1])-_b*x[idx+1];

  // differential eqn for h, the probability for the Na channel blocking
  // particle to be absent
  _a= 0.07*exp(-(x[idx]+37)/20);
  _b= 1/(exp(-0.1*(x[idx]+7))+1);
  dx[idx+2]= _a*(1.0-x[idx+2])-_b*x[idx+2];

  // differential eqn for n, the probability for one K channel activation
  // particle
  _a= -0.01*(x[idx]+27)/(exp(-0.1*(x[idx]+27))-1);
  _b= 0.125*exp(-(x[idx]+37)/80);
  dx[idx+3]= _a*(1.0-x[idx+3])-_b*x[idx+3];

  _a= 1/(0.15*(1+exp(-(x[idx]+38)/6.5)));
  _b= exp(-(x[idx]+38)/6.5)/(0.15*(1+exp(-(x[idx]+38)/6.5)));
  dx[idx+4]= _a*(1.0-x[idx+4])-_b*x[idx+4];
  
  // differential equation for the Ihf activation variable
  _a= 1/(1+exp((x[idx]+79.2)/9.78));
  _b= 0.51/(exp((x[idx]-1.7)/10)+exp(-(x[idx]+340)/52))+1;
  dx[idx+5]= (_a-x[idx+5])/_b;

  // differential equation for the Ihs activation variable
  _a= 1/(1+exp((x[idx]+71.3)/7.9));
  _b= 5.6/(exp((x[idx]-1.7)/14)+exp(-(x[idx]+260)/43))+1;
  dx[idx+6]= (_a-x[idx+6])/_b;
}

#endif
