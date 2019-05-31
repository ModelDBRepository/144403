/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#include "CN_demiGapsynapse.h"
#include "CN_synapse.cc"

demiGapsynapse:: demiGapsynapse(neuron* source, neuron *target, double ingsyn):
  synapse(source, target, demiGapIVARNO, demiGapPNO, DEMIGAPSYNAPSE)
{
  p[0]= ingsyn;               // injected current
}


demiGapsynapse::~demiGapsynapse()
{
}

double demiGapsynapse::gsyn()
{
  return p[0];
}

void demiGapsynapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double demiGapsynapse::Isyn(double *x)
{
  return p[0]*(source->E(x) - target->E(x));
}


