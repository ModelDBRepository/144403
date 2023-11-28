/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

// This is the base class for all types of plastic synapses with smooth
// changes in their maximal conductance. This class is abstract due to
// a lacking stdp function

#ifndef CN_ABSYNAPSE_SMSTDP_H
#define CN_ABSYNAPSE_SMSTDP_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define SMSTDPIVARNO 3

double SMSTDP_INIVARS[SMSTDPIVARNO]= {
  0.0,                    // 0 - S: internal variable
  1.0,                    // 1 - gsyn
  0.0                     // 2 - gsyn change rate
};

class absynapse_smSTDP: public synapse
{
 private:
  double tpre, tpost, dt;
  
 public:
  absynapse_smSTDP(neuron *, neuron *, int, int, int);
  ~absynapse_smSTDP();
  virtual double gsyn(double *);
  virtual double gsyn();
  virtual void set_gsyn(double, double *);
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void calc_dg(double *);
  virtual double stdp_fn(double) = 0;
  virtual void derivative(double *, double *);
};

#endif
