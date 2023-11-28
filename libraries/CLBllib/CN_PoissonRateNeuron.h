//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institute for Nonlinear Dynamics
//            University of California San Diego
//            La Jolla, CA 92093-0402
//
// email to:  tnowotny@ucsd.edu
//
// initial version: 2005-08-17
//
//--------------------------------------------------------------------------


#ifndef CN_POISSONRATENEURON_H
#define CN_POISSONRATENEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define POISSONRATENEURON_IVARNO 1

double *POISSONRATENEURON_p= PARAMS_ORN;

double POISSONRATENEURON_INIVARS[POISSONRATENEURON_IVARNO]= {
  -60.0                       // synspe variable
};

const char *POISSONRATENEURON_INIVARSTEXT[POISSONRATENEURON_IVARNO]= {
  "0 - membrane potential E"
};


// HH neuron class itself

class PoissonRateNeuron: public neuron
{
 private:

  double _a, _b;
  double tlast;


 public:
	  double Isyn;
  PoissonRateNeuron(int, double *);
  PoissonRateNeuron(int, vector<int>, double *);
  ~PoissonRateNeuron() { }
  inline virtual double E(double *);
  inline virtual double F(double *);
  virtual void validate_E(double *, double){}
  virtual void step(){}
  double S(double *);
  virtual void derivative(double *, double *);
  virtual void init(double *, double *);
  void SetVrest(double inVrest){mVrest = inVrest;}
  void ResetSynapse(double *x);
  double mVrest;
};

#endif



