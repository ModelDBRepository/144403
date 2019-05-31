/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
   Rallsynapse with inhibitory synapse learning rule of the Entorhnal cortex
   ... no saturation mechanism ... but delayed action of plasticity ...
--------------------------------------------------------------------------*/

#ifndef CN_RALLSYNAPSEECPLAST3_H
#define CN_RALLSYNAPSEECPLAST3_H

#include "queue.h"
#include "CN_Rallsynapse.h"

#define RALLECPLAST3PNO 6

#define RALLECPLAST3IVARNO 2

double RALLECPLAST3_INIVARS[RALLECPLAST3IVARNO]= {
  0.0,                    // 1 - S: internal variable
  0.0
};

class RallsynapseECplast3: public Rallsynapse
{  
 private:
  queue<double> chngTq;
  queue<double> dgq;
  double nextT;
  double nextdg;
  
 public:
  int synapse_change;
  RallsynapseECplast3(neuron *, neuron *, double, double, double, double,
			     double, double); 
  virtual ~RallsynapseECplast3();
  virtual void update_gsyn(double *);
  virtual double STDP_func(double);
};

#endif
