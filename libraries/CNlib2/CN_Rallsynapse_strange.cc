/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_RALLSYNAPSEST_CC
#define CN_RALLSYNAPSEST_CC

#include "CN_synapse.cc"
#define Heaviside(x) ((x > 0) ? 1 : 0)

// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag

RallsynapseST::RallsynapseST(neuron *insource, neuron *intarget,
			 double ingsyn, double inEsyn, double inEpre,
			 double intsyn, double inxmax,
			 int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= intsyn;           // synaptic timescale in msec
  p[4]= inxmax;           // maximal synaptic activation
} 

// This is the constructor to be used directly ...

RallsynapseST::RallsynapseST(neuron *insource, neuron *intarget,
			 double ingsyn, double inEsyn, double inEpre,
			 double intsyn, double inxmax):
  synapse(insource, intarget, RIVARNOST, RPNOST, RALL)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= intsyn;           // synaptic timescale in msec
  p[4]= inxmax;           // maximal synaptic activation
} 

RallsynapseST::RallsynapseST(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, RIVARNOST, RPNOST, RALL)
{
  set_p(inp);
} 

RallsynapseST::~RallsynapseST()
{
}

double RallsynapseST::gsyn()
{
  return p[0];
}

void RallsynapseST::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double RallsynapseST::Isyn(double *x)
{
  return -p[0]*(p[4]-x[idx+1])*(target->E(x)-p[1]);
}


void RallsynapseST::derivative(double *x, double *dx)
{
  dx[idx+0]= (-x[idx]+Heaviside(source->E(x)-p[2]))/p[3];
  dx[idx+1]= ((p[4]-x[idx+1])/2.0-x[idx])*x[idx+1]/(p[3]*p[4]);
}

// end of class implementation

#endif


