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


#ifndef CN_VDPOLNEURON_H
#define CN_VDPOLNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define VDPOL_IVARNO 2
#define VDPOL_PNO 3

double stdVDPOL_p[VDPOL_PNO]= {
  1.0,          // 0 - eta
  0.1,          // 1 - omega^2
  0.0           // 2 - noise level
};

double *VDPOL_p= stdVDPOL_p;

const char *VDPOL_p_text[VDPOL_PNO]= {
  "0 - eta",
  "1 - omega^2",
  "2 - noise level"
};

double VDPOL_INIVARS[VDPOL_IVARNO]= {
  0.1,                       // 0 - amplitude
  0.0                        // 1 - internal var
};

const char *VDPOL_INIVARSTEXT[VDPOL_IVARNO]= {
  "0 - amplitude",
  "1 - internal var"
};


// VDPOL neuron class itself

class VdPolneuron: public neuron
{
 private:
  double Isyn;
 public:
  VdPolneuron(int, double *);
  VdPolneuron(int, vector<int>, double *);
  ~VdPolneuron() { }
  inline virtual double E(double *);
  virtual void derivative(double *, double *);
  virtual void noise(double *, double *);
};

#endif



