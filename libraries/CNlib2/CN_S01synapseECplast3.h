/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
   S01synapse with inhibitory synapse learning rule of the Entorhinal cortex
   ... no saturation mechanism ... but delayed action of plasticity ...
--------------------------------------------------------------------------*/

#ifndef CN_S01SYNAPSEECPLAST3_H
#define CN_S01SYNAPSEECPLAST3_H

#include "queue.h"
#include "CN_S01synapse.h"

#define S01ECPLAST3PNO 8

#define S01ECPLAST3IVARNO 1

double S01ECPLAST3_INIVARS[S01ECPLAST3IVARNO]= {
  0.0                    // 1 - S: internal variable
};

class S01synapseECplast3: public S01synapse
{  
 private:
  queue<double> chngTq;
  queue<double> dgq;
  double nextT;
  double nextdg;
  
 public:
  int synapse_change;
  S01synapseECplast3(neuron *, neuron *, double, double, double, double,
			     double, double, double, double); 
  virtual ~S01synapseECplast3();
  virtual void update_gsyn(double *);
  virtual double STDP_func(double);
};

#endif
