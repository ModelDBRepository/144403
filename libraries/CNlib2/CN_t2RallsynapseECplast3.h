/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
   T2Rallsynapse with inhibitory synapse learning rule of the Entorhnal cortex
   ... no saturation mechanism ... but delayed action of plasticity ...
--------------------------------------------------------------------------*/

#ifndef CN_T2RALLSYNAPSEECPLAST3_H
#define CN_T2RALLSYNAPSEECPLAST3_H

#include "queue.h"
#include "CN_t2Rallsynapse.h"

#define T2RALLECPLAST3PNO 7

#define T2RALLECPLAST3IVARNO 2

double T2RALLECPLAST3_INIVARS[T2RALLECPLAST3IVARNO]= {
  0.0,                    // 1 - S: internal variable
  0.0
};

class t2RallsynapseECplast3: public t2Rallsynapse
{  
 private:
  queue<double> chngTq;
  queue<double> dgq;
  double nextT;
  double nextdg;
  
 public:
  int synapse_change;
  t2RallsynapseECplast3(neuron *, neuron *, double, double, double, double,
			     double, double, double); 
  virtual ~t2RallsynapseECplast3();
  virtual void update_gsyn(double *);
  virtual double STDP_func(double);
};

#endif
