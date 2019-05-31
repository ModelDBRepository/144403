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


#ifndef CN_DATA_H
#define CN_DATA_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define Data_IVARNO 0
#define Data_PNO 2

double stdData_p[Data_PNO]= {
  0.0, // t0
  0.1 // dt
};

double *Data_p= stdData_p;

const char *Data_p_text[Data_PNO]= {
  "t0",
  "dt"
};



class Data: public neuron
{
 private:
  int cnt;
  double *dA;
 public:
  Data(int, char *, double *);
  virtual Data::~Data();
  inline virtual double E(double *);
  virtual void derivative(double *, double *) { }
};

#endif



