/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_T2RALLSYNAPSE_CC
#define CN_T2RALLSYNAPSE_CC

#include "CN_synapse.cc"
#define Heaviside(x) ((x > 0) ? 1 : 0)

// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag

t2Rallsynapse::t2Rallsynapse(neuron *insource, neuron *intarget,
			 double ingsyn, double inEsyn, double inEpre,
			 double intauR, double intauF,
			 int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= intauR;           // synaptic timescale in msec
  p[4]= intauF;
} 

// This is the constructor to be used directly ...

t2Rallsynapse::t2Rallsynapse(neuron *insource, neuron *intarget,
			 double ingsyn, double inEsyn, double inEpre,
			 double intauR, double intauF):
  synapse(insource, intarget, T2RIVARNO, T2RPNO, T2RALL)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= intauR;           // synaptic timescale in msec
  p[4]= intauF;
} 

t2Rallsynapse::t2Rallsynapse(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, T2RIVARNO, T2RPNO, T2RALL)
{
  set_p(inp);
} 

t2Rallsynapse::~t2Rallsynapse()
{
}

double t2Rallsynapse::gsyn()
{
  return p[0];
}

void t2Rallsynapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double t2Rallsynapse::Isyn(double *x)
{
  return -p[0]*x[idx+1]*(target->E(x)-p[1]);
}


void t2Rallsynapse::derivative(double *x, double *dx)
{
  dx[idx+0]= (-x[idx]+Heaviside(source->E(x)-p[2]))/p[3];
  dx[idx+1]= -x[idx+1]/p[4]+x[idx]/p[3];
}

// end of class implementation

#endif


