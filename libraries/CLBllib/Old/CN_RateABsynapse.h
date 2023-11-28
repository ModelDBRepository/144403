/*--------------------------------------------------------------------------
   Author: Thomas Nowotny

   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402

   email to:  tnowotny@ucsd.edu

   initial version: 2005-08-17

--------------------------------------------------------------------------*/


#ifndef CN_RATEABSYNAPSE_H
#define CN_RATEABSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define RATEABSYNPNO 7
#define RATEABSYNIVARNO 1

double RATEABSYN_INIVARS[RATEABSYNIVARNO]= {
  0.0                    // 0 - S: internal variable
};

class RateABsynapse: public synapse
{
 private:

 public:
  RateABsynapse(neuron *, neuron *, double, double, double, double, double,
	    double, double, int, int, int);
  RateABsynapse(neuron *, neuron *, double, double, double, double, double,
	    double, double);
  RateABsynapse(neuron *, neuron *, double *);
  virtual ~RateABsynapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void init(double *, double *);
  virtual void derivative(double *, double *);
  virtual void noise(double *, double *);
  double S(double *);
};

#endif
