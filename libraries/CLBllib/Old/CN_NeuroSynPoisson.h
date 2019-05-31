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


#ifndef CN_NEUROSYNPOISSON_H
#define CN_NEUROSYNPOISSON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

double *NEUROSYNPOISSON_p= PARAMS_LN;


// HH neuron class itself

class NeuroSynPoisson: public neuron
{
 private:
  double Isyn;
  double ICa;
  double IKCa;
  double _a, _b;
  double tlast;
  int refract;
  double spikingThistime;
 public:
  NeuroSynPoisson(int, double *);
  NeuroSynPoisson(int, vector<int>, double *);
  ~NeuroSynPoisson() { }
  inline virtual double E(double *);
  inline virtual double F(double *){return 0;}
  double S(double *);
  virtual void derivative(double *, double *);
  virtual void init(double *, double *);
  double Transfer(double isyn){return 0.0;}
  void ResetSynapse(double *x);
  void advance(double *, double,double);
};

#endif



