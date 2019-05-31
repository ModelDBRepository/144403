/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_S01SYNAPSE_CC
#define CN_S01SYNAPSE_CC

#include "CN_synapse.cc"

// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag

S01synapse::S01synapse(neuron *insource, neuron *intarget,
		     double ingsyn, double inEsyn, double inEpre,
		     double intau, double inS1, double inVslope,
		     int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= intau;           // alpha timescale in 1/msec
  p[4]= inS1;           // beta timescale in 1/msec
  p[5]= inVslope;         // steepness of activation curve as func of Vpre  
} 

// This is the constructor to be used directly ...

S01synapse::S01synapse(neuron *insource, neuron *intarget,
		     double ingsyn, double inEsyn, double inEpre,
		     double intau, double inS1, double inVslope):
  synapse(insource, intarget, S01SYNIVARNO, S01SYNPNO, S01SYN)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= intau;           // timescale in 1/msec
  p[4]= inS1;           // S_1 
  p[5]= inVslope;         // steepness of activation curve as func of Vpre  
} 

S01synapse::S01synapse(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, S01SYNIVARNO, S01SYNPNO, S01SYN)
{
  set_p(inp);
} 

S01synapse::~S01synapse()
{
}

double S01synapse::gsyn()
{
  return p[0];
}

void S01synapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double S01synapse::Isyn(double *x)
{
  return -p[0]*x[idx]*(target->E(x)-p[1]);
}


void S01synapse::derivative(double *x, double *dx)
{
  static double tmp;
  
  tmp= (1.0+tanh((source->E(x)-p[2])/p[5]))/2.0;
  dx[idx]= (tmp - x[idx])/(p[3]*(p[4]-tmp));
}

// end of class implementation

#endif


