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
#define CN_ABSYNAPSE_SMSTDP1_H

#include "CN_absynapse_smSTDP.h"
#include "CN_neuron.h"

#define SMSTDP1PNO 13

double stdSMSTDP1_p[SMSTDP1PNO]= {
  0.0,   // 0 - Erev: reversal potential
  -20.0, // 1 - Esyn: presynaptic threshold
  0.5,   // 2 - asyn: alpha (rise rate)
  0.1,   // 3 - bsyn: beta (fall rate)
  10.0,  // 4 - Vslope: Slope of activation curve
  0.01,  // 5 - decay: rate of STDP action time window
  10.0,  // 6 - gmax: maximal g
  0.01,  // 7 - g0: equilibrium g
  0.001, // 8 - gdecay: time scale of g decay towards g0
  0.01,   // 9 - Aplus: amplitude of positive STDP curve
  0.01,   // 10 - Aminus: amplitude of negative STDP curve
  20.0,  // 11 - tauplus: time scale of positive STDP curve
  30.0   // 12 - tauminus: time scale of negative STDP curve
};

class absynapse_smSTDP1: public absynapse_smSTDP
{
 public:
  absynapse_smSTDP1(neuron *, neuron *, double, double, double,
		    double, double, double, double, double, double, double,
		    double, double, double);
  absynapse_smSTDP1(neuron *, neuron *, double *);
  ~absynapse_smSTDP1();
  virtual double stdp_fn(double);
};

#endif
