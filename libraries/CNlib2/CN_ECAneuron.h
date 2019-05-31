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


#ifndef DC_ECNEURON_H
#define DC_ECNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the Entorhinal cortex stellate cell

#define ECN_IVARNO 7
#define ECN_PNO 10

double stdECN_p[ECN_PNO]= {
  52.0,          // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
  55.0,          // 1 - ENa: Na equi potential in mV
  0.5,           // 2 - gnap
  11.0,          // 3 - gK: K conductance in 1/(mOhms * cm^2)
  -90.0,         // 4 - EK: K equi potential in mV
  0.5,           // 5 - gl: leak conductance in 1/(mOhms * cm^2)
  -65.0,         // 6 - El: leak equi potential in mV
  1.5,           // 7 - gh
  -20.0,         // 8 - Vh
  1.5            // 9 - Cmem: membr. capacity density in muF/cm^2
};

double *ECN_p= stdECN_p;

const char *ECN_p_text[ECN_PNO]= {
  "0 - gNa: Na conductance in 1/(mOhms * cm^2)",
  "1 - ENa: Na equi potential in mV",
  "2 - gNap",
  "3 - gK: K conductance in 1/(mOhms * cm^2)",
  "4 - EK: K equi potential in mV",
  "5 - gl: leak conductance in 1/(mOhms * cm^2)",
  "6 - El: leak equi potential in mV",
  "7 - gh",
  "8 - Vh",
  "9 - Cmem: membr. capacity density in muF/cm^2"
};

double ECN_INIVARS[ECN_IVARNO]= {
  -53.77902178,                   // 0 - membrane potential E
  0.0262406368,                   // 1 - prob. for Na channel activation m
  0.9461831106,                   // 2 - prob. for not Na channel blocking h
  0.1135915933,                   // 3 - prob. for K channel activation n
  0.08109646237,                  // 4 - Nap
  0.06918464221,                  // 5 - Ih1 activation
  0.09815937825                   // 6 - Ih2 activation  
};

const char *ECN_INIVARSTEXT[ECN_IVARNO]= {
  "0 - membrane potential E",
  "1 - prob. for Na channel activation m",
  "2 - prob. for not Na channel blocking h",
  "3 - prob. for K channel activation n",
  "4 - m_Nap",
  "5 - Ih1 activation",
  "6 - Ih2 activation"
};


// ECneuron class itself

class ECneuron: public neuron
{
 private:
  double Isyn;
  double _a, _b;
 public:
  ECneuron(int, double *);
  ECneuron(int, vector<int>, double *);
  ~ECneuron() { }
  inline virtual double E(double *);
  virtual void derivative(double *, double *);
};

#endif



