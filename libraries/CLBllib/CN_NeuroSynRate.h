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


#ifndef CN_NEUROSYNRATE_H
#define CN_NEUROSYNRATE_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define NEUROSYNRATE_IVARNO 1

double *NEUROSYNRATE_p= PARAMS_LN;

double NEUROSYNRATE_INIVARS[NEUROSYNRATE_IVARNO]= {
  -60.0                       // synspe variable
};

const char *NEUROSYNRATE_INIVARSTEXT[NEUROSYNRATE_IVARNO]= {
  "0 - membrane potential E"
};


// HH neuron class itself

class NeuroSynRate: public neuron
{
 private:

  double ICa;
  double IKCa;
  double _a, _b;
  double tlast;
  double mF;


 public:
	  double Isyn;
  NeuroSynRate(int, double *);
  NeuroSynRate(int, vector<int>, double *);
  ~NeuroSynRate() { }
  inline virtual double E(double *);
  inline virtual double F(double *);
  double S(double *);
  virtual void derivative(double *, double *);
  virtual void init(double *, double *);
  double Transfer(double isyn);
  void SetVrest(double inVrest){mVrest = inVrest;}
  void ResetSynapse(double *x);
  double mVrest;
};

#endif



