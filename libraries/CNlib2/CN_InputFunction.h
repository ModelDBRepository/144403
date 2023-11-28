/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_INPUTFUNCTION_H
#define CN_INPUTFUNCTION_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define DCIVARNO 0
#define DCPNO 2

class InputFunction: public synapse
{
 public:
  InputFunction(neuron *, double, double);
  virtual ~InputFunction();
  virtual void set_A(double);
  virtual void set_Om(double);
  virtual double gsyn() {return 1.0;}
  virtual void set_gsyn(double) { }
  virtual double Isyn(double *);
  virtual void derivative(double *, double *) { }
};

#endif

