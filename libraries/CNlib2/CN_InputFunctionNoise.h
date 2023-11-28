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
#define DCPNO 3

class InputFunctionNoise: public synapse
{
 public:
  InputFunctionNoise(neuron *, double, double, double);
  virtual ~InputFunctionNoise();
  virtual void set_A(double);
  virtual void set_Om(double);
  virtual void set_Noise(double);
  virtual double gsyn() {return 1.0;}
  virtual void set_gsyn(double) { }
  virtual double Isyn(double *);
  virtual void derivative(double *, double *) { }
};

#endif

