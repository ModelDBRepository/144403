/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_ECDEMIGAPSYNAPSE_H
#define CN_ECDEMIGAPSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define ECdemiGapPNO 4
#define ECdemiGapIVARNO 1

double ECdemiGap_INIVARS[ECdemiGapIVARNO]= {
  0.25                    // 0 - S: internal variable
};


class ECdemiGapsynapse: public synapse
{
 public:
  ECdemiGapsynapse(neuron *, neuron *, double, double, double, double);
  virtual ~ECdemiGapsynapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *);
};

#endif
