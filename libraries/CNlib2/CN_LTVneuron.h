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


#ifndef CN_LTVNEURON_H
#define CN_LTVNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define LTV_IVARNO 1
#define LTV_PNO 2

double stdLTV_p[LTV_PNO]= {
  1.0,          // 0 - rho_ii: "self inhibition"
  0.0           // 1 - sigma: noise level
};

double *LTV_p= stdLTV_p;

const char *LTV_p_text[LTV_PNO]= {
  "0 - rho_ii: self inhibition",
};

double LTV_INIVARS[LTV_IVARNO]= {
  0.1                       // 0 - firing rate
};

const char *LTV_INIVARSTEXT[LTV_IVARNO]= {
  "0 - firing rate"
};


// LTV neuron class itself

class LTVneuron: public neuron
{
 private:
  double Isyn;
 public:
  LTVneuron(int, double *);
  LTVneuron(int, vector<int>, double *);
  ~LTVneuron() { }
  inline virtual double E(double *);
  virtual void derivative(double *, double *);
  virtual void noise(double *, double *);
};

#endif



