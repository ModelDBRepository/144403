/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_RALLSYNAPSEST_H
#define CN_RALLSYNAPSEST_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define RPNOST 5

#define RIVARNOST 2

double RASYNST_INIVARS[RIVARNOST]= {
  0.0,                    // 0 - f: internal variable 
  0.0                     // 1 - g: internal variable
};

class RallsynapseST: public synapse
{
  
 public:
  RallsynapseST(neuron *, neuron *, double, double, double, double, double,
	      int, int, int);
  RallsynapseST(neuron *, neuron *, double, double, double, double, double);
  RallsynapseST(neuron *, neuron *, double *);
  virtual ~RallsynapseST();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *);
};

#endif
