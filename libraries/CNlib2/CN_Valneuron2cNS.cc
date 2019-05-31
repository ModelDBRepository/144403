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

Valneuron::Valneuron(int inlabel, double *the_p= Val_p):
  neuron(inlabel, Val_IVARNO, VALNEURON, the_p, Val_PNO)
{
}

Valneuron::Valneuron(int inlabel, vector<int> inpos, double *the_p= Val_p):
  neuron(inlabel, Val_IVARNO, VALNEURON, inpos, the_p, Val_PNO)
{
}

// note that this is the dendritic memb. potential!
inline double Valneuron::E(double *x)
{
  assert(enabled);
  return x[idx+4];
}

// soma
inline double Valneuron::Esoma(double *x)
{
  assert(enabled);
  return x[idx];
}


void Valneuron::derivative(double *x, double *dx)
{
  static double IVV;

  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }
  
  IVV= p[8]*(x[idx]-x[idx+4]);
  // differential eqn for E, the membrane potential
  dx[idx]= -(pw3(x[idx+1])*x[idx+2]*p[0]*(x[idx]-p[1]) +
	     pw4(x[idx+3])*p[2]*(x[idx]-p[3])+
	     p[4]*(x[idx]-p[5])+IVV-p[11])/p[6];
  
  // diferential eqn for m, the probability for one Na channel activation
  // particle
  _a= 0.32*(-52.0-x[idx]) / (exp((-52.0-x[idx])/4.0)-1.0);
  _b= 0.28*(x[idx]+25)/(exp((x[idx]+25)/5.0)-1.0);
  dx[idx+1]= _a*(1.0-x[idx+1])-_b*x[idx+1];

  // differential eqn for h, the probability for the Na channel blocking
  // particle to be absent
  _a= 0.128*exp((-48-x[idx])/18.0);   
  _b= 4.0 / (exp((-25-x[idx])/5.0)+1.0);
  dx[idx+2]= _a*(1.0-x[idx+2])-_b*x[idx+2];

  // differential eqn for n, the probability for one K channel activation
  // particle
  _a= .032*(-50-x[idx]) / (exp((-50.0-x[idx])/5.0)-1.0); 
  _b= 0.5*exp((-55-x[idx])/40.0);
  dx[idx+3]= _a*(1.0-x[idx+3])-_b*x[idx+3];

  dx[idx+4]= (Isyn + IVV - p[10]*(x[idx+4]-p[5]))/p[9];
}

void Valneuron::noise(double *x, double *dx)
{
  dx[idx]= RG.n()*p[10];
  for (int i= 1; i < iVarNo; i++) {
    dx[idx+i]= 0.0;
  }
}

// overloading this one to work on Esoma (!)
void Valneuron::spike_detect(double *x)
{
  assert(enabled);
  if (Esoma(x) >= SPK_V_THRESH)
  {
    if (!spiking)
    {
      start_spiking= 1;
      spiking= 1;
      spike_time= x[0];
    }
    else start_spiking= 0;
  }
  else {
    spiking= 0;
    start_spiking= 0;
  }
}

#endif
