/*--------------------------------------------------------------------------
   Author: Thomas Nowotny

   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402

   email to:  tnowotny@ucsd.edu

   initial version: 2005-08-17

--------------------------------------------------------------------------*/

#ifndef CN_NEUROSYN_CC
#define CN_NEUROSYN_CC

#include "CN_neuron.cc"

NeuroSyn::NeuroSyn(int inlabel, double *the_p= NEUROSYN_p):
  neuron(inlabel, NEUROSYN_IVARNO, NEUROSYN, the_p, NEUROSYN_PNO)
{
	 tlast= -10000;
}

NeuroSyn::NeuroSyn(int inlabel, vector<int> inpos, double *the_p= NEUROSYN_p):
  neuron(inlabel, NEUROSYN_IVARNO, NEUROSYN, inpos, the_p, NEUROSYN_PNO)
{
	 tlast= -10000;
}

inline double NeuroSyn::E(double *x)
{
  assert(enabled);
  return x[idx];
}

double NeuroSyn::S(double *x)
{
  assert(enabled);
  return x[idx+5];
}
void NeuroSyn::derivative(double *x, double *dx)
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

  static double dt;

   dt= x[0] - tlast;
   if ((dt >= 0) && (dt <= p[16])) {
 	  dx[idx+5]= p[14] - p[15]*x[idx+5];
   } else {
     if ((x[idx] > p[13]) && (dt > p[16])) {
       // new spike ... start releasing
       tlast= x[0];
      dx[idx+5]= p[14] - p[15]*x[idx+5];
     }
     else {
       // no release
       dx[idx+5]= -p[15]*x[idx+5];
     }
   }
}

void NeuroSyn::init(double *x, double *iniVars)
{

	assert(enabled);
	  for (int i= 0; i < iVarNo; i++)
	    x[idx+i]= iniVars[i];

	  start_spiking= 0;
	  spiking= 0;
	  spike_time= -1.0;

	  tlast= -10000;
}

void NeuroSyn::ResetSynapse(double *x)
{
	x[idx+5] = 0;
}

#endif
