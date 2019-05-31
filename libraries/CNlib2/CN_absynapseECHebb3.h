/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
   absynapse with inhibitory synapse learning rule of the Entorhnal cortex
   ... no saturation mechanism ... but delayed action of plasticity ...
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSEECHEBB3_H
#define CN_ABSYNAPSEECHEBB3_H

#include "queue.h"
#include "CN_absynapse.h"

#define ABECHEBB3PNO 8

#define ABECHEBB3IVARNO 1

double ABECHEBB3_INIVARS[ABECHEBB3IVARNO]= {
  0.0                    // 1 - S: internal variable
};

class absynapseECHebb3: public absynapse
{  
 private:
  queue<double> chngTq;
  queue<double> dgq;
  double nextT;
  double nextdg;
  
 public:
  int synapse_change;
  absynapseECHebb3(neuron *, neuron *, double, double, double, double,
			     double, double, double, double); 
  virtual ~absynapseECHebb3();
  virtual void update_gsyn(double *);
  virtual double STDP_func(double);
};

#endif
