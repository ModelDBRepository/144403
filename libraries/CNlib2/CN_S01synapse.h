/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_S01SYNAPSE_H
#define CN_S01SYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define S01SYNPNO 6
#define S01SYNIVARNO 1

double S01SYN_INIVARS[S01SYNIVARNO]= {
  0.0                    // 0 - S: internal variable 
};

class S01synapse: public synapse
{
  
 public:
  S01synapse(neuron *, neuron *, double, double, double, double, double,
	    double, int, int, int);
  S01synapse(neuron *, neuron *, double, double, double, double, double,
	    double);
  S01synapse(neuron *, neuron *, double *);
  virtual ~S01synapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *);
};

#endif
