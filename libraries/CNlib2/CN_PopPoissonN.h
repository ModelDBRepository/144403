/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_POPPOISSONN_H
#define CN_POPPOISSONN_H

#include "CN_neuron.h"

#define POPPOI_IVARNO 0
#define POPPOI_PNO 4

double stdPOPPOI_p[POPPOI_PNO]= {
  0.0,  // refractory period [ms]
  -65.0, // resting potential
  20.0, // spiking potential
  0.5 // spike duration
};

double *POPPOI_p= stdPOPPOI_p;

class PopPoissonN: public neuron
{
 private:
  double fire_t;
  double next_spike;
  double *Lambda;
  
 public:
  PopPoissonN(int, vector<int>, double *);
  PopPoissonN(int, double *);
  ~PopPoissonN();
  virtual double E(double *);
  virtual void validate_E(double *, double);
  virtual void step();
  virtual void derivative(double *, double *) { }
  virtual void init(double *, double *);
  virtual void set_Lambda(double *);  // set an external rate variable
};

#endif





