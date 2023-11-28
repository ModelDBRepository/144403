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


#ifndef CN_NEUROSYNADAPT_H
#define CN_NEUROSYNADAPT_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

double *NEUROSYN_p= PARAMS_LN;


// HH neuron class itself

class NeuroSynAdapt: public neuron
{
 private:
  double ICa;
  double IKCa;
  double _a, _b;
  double tlast;
 public:
	  double Isyn;
  NeuroSynAdapt(int, double *);
  NeuroSynAdapt(int, vector<int>, double *);
  ~NeuroSynAdapt() { }
  inline virtual double E(double *);
  inline virtual double F(double *){return 0;}
  double S(double *);
  double Ca(double *);
  double Theta(double *);
  virtual void derivative(double *, double *);
  virtual void init(double *, double *);
  double Transfer(double isyn){return 0.0;}
  void ResetSynapse(double *x);
};

#endif



