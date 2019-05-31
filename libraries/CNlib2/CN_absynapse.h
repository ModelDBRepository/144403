/*--------------------------------------------------------------------------
   Author: Thomas Nowotny

   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402

   email to:  tnowotny@ucsd.edu

   initial version: 2005-08-17

--------------------------------------------------------------------------*/


#ifndef CN_ABSYNAPSE_H
#define CN_ABSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define ABSYNPNO 7
#define ABSYNIVARNO 1

double ABSYN_INIVARS[ABSYNIVARNO]= {
  0.0                    // 0 - S: internal variable
};

class absynapse: public synapse
{
 private:
  double tlast;

 public:
  absynapse(neuron *, neuron *, double, double, double, double, double,
	    double, double, int, int, int);
  absynapse(neuron *, neuron *, double, double, double, double, double,
	    double, double);
  absynapse(neuron *, neuron *, double *);
  virtual ~absynapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void init(double *, double *);
  virtual void derivative(double *, double *);
  virtual void noise(double *, double *);
  double S(double *);
};

#endif
