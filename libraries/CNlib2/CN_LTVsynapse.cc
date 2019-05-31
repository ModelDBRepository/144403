/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_LTVSYNAPSE_CC
#define CN_LTVSYNAPSE_CC

#include "CN_synapse.cc"

// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag

LTVsynapse::LTVsynapse(neuron *insource, neuron *intarget, double ingsyn,
		       int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
} 

// This is the constructor to be used directly ...

LTVsynapse::LTVsynapse(neuron *insource, neuron *intarget, double ingsyn):
  synapse(insource, intarget, LTVIVARNO, LTVPNO, LTVSYN)
{
  p[0]= ingsyn;           // gsyn strength of synapse
} 

LTVsynapse::LTVsynapse(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, LTVIVARNO, LTVPNO, LTVSYN)
{
  set_p(inp);
} 

LTVsynapse::~LTVsynapse()
{
}

double LTVsynapse::gsyn()
{
  return p[0];
}

void LTVsynapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double LTVsynapse::Isyn(double *x)
{
  return p[0]*source->E(x)*target->E(x);
}

// end of class implementation

#endif


