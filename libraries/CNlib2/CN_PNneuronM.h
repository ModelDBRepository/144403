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


#ifndef CN_PNNEURON_H
#define CN_PNNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define PN_IVARNO 5
#define PN_PNO 12

double stdPN_p[PN_PNO]= {
		7.15, // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
			50.0, // 1 - ENa: Na equi potential in mV
			1.43, // 2 - gK: K conductance in 1/(mOhms * cm^2)
			-95.0, // 3 - EK: K equi potential in mV
			0.021, // 4 - gl: leak conductance in 1/(mOhms * cm^2)
			-55.0, // 5 - El: leak equi potential in mV
			0.00572, // 6 - gKl: potassium leakage conductivity
			-95.0, // 7 - EKl: potassium leakage equi pot in mV
			65, // 8 - V0: ~ total equi potential (?)
			0.143, // 9 - Cmem: membr. capacity density in muF/cm^2
			0.715, // 10 - gM: conductance of the M current
			0.041270138525220 // 11- IDC: baseline offset current
};




double *PN_p= stdPN_p;

const char *PN_p_text[PN_PNO]= {
  "0 - gNa: Na conductance in 1/(mOhms * cm^2)",
  "1 - ENa: Na equi potential in mV",
  "2 - gK: K conductance in 1/(mOhms * cm^2)",
  "3 - EK: K equi potential in mV",
  "4 - gl: leak conductance in 1/(mOhms * cm^2)",
  "5 - El: leak equi potential in mV",
  "6 - gKl: potassium leakage conductivity",
  "7 - EKl: potassium leakage equi pot in mV",
  "8 - V0: ~ total equi potential (?)",
  "9 - Cmem: membr. capacity density in muF/cm^2",
  "10 - gM: conductance of the M current",
  "11 - IDC: baseline offset current"
};

double PN_INIVARS[PN_IVARNO]= {
  -60.0,                       // 0 - membrane potential E
  0.0529324,                   // 1 - prob. for Na channel activation m
  0.3176767,                   // 2 - prob. for not Na channel blocking h
  0.5961207,                   // 3 - prob. for K channel activation n
  0.1                          // 4 - M current activation
};

const char *PN_INIVARSTEXT[PN_IVARNO]= {
  "0 - membrane potential E",
  "1 - prob. for Na channel activation m",
  "2 - prob. for not Na channel blocking h",
  "3 - prob. for K channel activation n",
  "4 - M current activation"
};


// HH neuron class itself

class PNneuron: public neuron
{
 private:
  double Isyn;
  double ICa;
  double IKCa;
  double _a, _b;
 public:
  PNneuron(int, double *);
  PNneuron(int, vector<int>, double *);
  ~PNneuron() { }
  inline virtual double E(double *);
  inline virtual double F(double *){return 0;}
  inline virtual double S(double *){return 0;}
  virtual void derivative(double *, double *);
};

#endif



