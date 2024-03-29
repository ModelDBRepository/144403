/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_DEMIGAPSYNAPSE_H
#define CN_DEMIGAPSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define demiGapPNO 1
#define demiGapIVARNO 0

class demiGapsynapse: public synapse
{
 public:
  demiGapsynapse(neuron *, neuron *, double);
  virtual ~demiGapsynapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *) { }
};

#endif
