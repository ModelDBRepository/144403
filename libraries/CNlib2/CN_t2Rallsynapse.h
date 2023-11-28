/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_T2RALLSYNAPSE_H
#define CN_T2RALLSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define T2RPNO 5

#define T2RIVARNO 2

double RASYN_INIVARS[T2RIVARNO]= {
  0.0,                    // 0 - f: internal variable 
  0.0                     // 1 - g: internal variable
};

class t2Rallsynapse: public synapse
{
  
 public:
  t2Rallsynapse(neuron *, neuron *, double, double, double, double, double,
	      int, int, int);
  t2Rallsynapse(neuron *, neuron *, double, double, double, double, double);
  t2Rallsynapse(neuron *, neuron *, double *);
  virtual ~t2Rallsynapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *);
};

#endif
