/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_LTVSYNAPSE_H
#define CN_LTVSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define LTVPNO 1

#define LTVIVARNO 0

class LTVsynapse: public synapse
{
  
 public:
  LTVsynapse(neuron *, neuron *, double, int, int, int);
  LTVsynapse(neuron *, neuron *, double);
  LTVsynapse(neuron *, neuron *, double *);
  virtual ~LTVsynapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *) { }
};

#endif
