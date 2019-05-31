/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSE_SMSTDP_CC
#define CN_ABSYNAPSE_SMSTDP_CC

#include "CN_synapse.cc"

// There is no constructor to be used directly (class is abstract) ...

absynapse_smSTDP::absynapse_smSTDP(neuron *insource, neuron *intarget,
				     int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  tpre= -1e10;
  tpost= -1e10;
  dt= 0.0;
} 

absynapse_smSTDP::~absynapse_smSTDP()
{
}

double absynapse_smSTDP::gsyn(double *x)
{
  return x[idx+1];
}

// this function does not work ...
double absynapse_smSTDP::gsyn()
{
  cerr << "using unsupported function!" << endl;
  exit(1);
}

void absynapse_smSTDP::set_gsyn(double ingsyn, double *x)
{
  x[idx+1]= ingsyn;
}

void absynapse_smSTDP::set_gsyn(double ingsyn)
{
  cerr << "using unsupported function!" << endl;
  exit(1);
}

double absynapse_smSTDP::Isyn(double *x)
{
  return -x[idx+1]*x[idx]*(target->E(x)-p[0]);
}

void absynapse_smSTDP::calc_dg(double *x)
{
  if (source->start_spiking) {
    tpre= x[0];
    dt= tpost-tpre;
    x[idx+2]+= stdp_fn(dt);
  }
  if (target->start_spiking) {
    tpost= x[0];
    dt= tpost-tpre;
    x[idx+2]+= stdp_fn(dt);
  }
}

void absynapse_smSTDP::derivative(double *x, double *dx)
{
  dx[idx]= p[2]*(1.0-x[idx])*(1.0+tanh((source->E(x)-p[1])/p[4]))/2.0
    -p[3]*x[idx];
  dx[idx+1]= x[idx+2]*x[idx+1]*(p[6]-x[idx+1])-p[8]*(x[idx+1]-p[7]);
  dx[idx+2]= -p[5]*x[idx+2];
}

// end of class implementation

#endif


