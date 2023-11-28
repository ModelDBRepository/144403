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


#ifndef CN_COLPITTS_H
#define CN_COLPITTS_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define Colp_IVARNO 3
#define Colp_PNO 4

double stdColp_p[Colp_PNO]= {
  1.0, // a
  0.0797, // g
  0.6898, // q
  6.2723 // eta
};

double *Colp_p= stdColp_p;
  
const char *Colp_p_text[Colp_PNO]= {
  "a",
  "g",
  "q",
  "eta"
};

double Colp_INIVARS[Colp_IVARNO]= {
  0.02,
  0.69,
  -0.53
};

const char *Colp_INIVARSTEXT[Colp_IVARNO]= {
  "x0",
  "x1",
  "x2"
};


// Colpitts oscillator

class Colpitts: public neuron
{
 private:
  double Isyn;

 public:
  Colpitts(int, double *);
  Colpitts(int, vector<int>, double *);
  ~Colpitts() { }
  inline virtual double E(double *);
  virtual void derivative(double *, double *);
};

#endif



