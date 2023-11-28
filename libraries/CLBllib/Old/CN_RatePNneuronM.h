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


#ifndef CN_RATEPNNEURON_H
#define CN_RATEPNNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define RATEPN_IVARNO 0
#define RATEPN_PNO 12

double stdRATEPN_p[RATEPN_PNO] = {0, // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
		0, // 1 - ENa: Na equi potential in mV
		0, // 2 - gK: K conductance in 1/(mOhms * cm^2)
		0, // 3 - EK: K equi potential in mV
		0, // 4 - gl: leak conductance in 1/(mOhms * cm^2)
		0, // 5 - El: leak equi potential in mV
		0., // 6 - gKl: potassium leakage conductivity
        0, // 7 - EKl: potassium leakage equi pot in mV
		0, // 8 - V0: ~ total equi potential (?)
		0.143, // 9 - Cmem: membr. capacity density in muF/cm^2
		0.715, // 10 - gM: conductance of the M current
		0.041270138525220 // 11- IDC: baseline offset current
		};

double *RATEPN_p= stdRATEPN_p;

const char *RATEPN_p_text[RATEPN_PNO]= {
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




// HH neuron class itself

class RatePNneuron: public neuron
{
 private:
  double Isyn;
  double ICa;
  double IKCa;
  double _a, _b;
 public:
	 RatePNneuron(int, double *);
	 RatePNneuron(int, vector<int>, double *);
  ~ RatePNneuron() { }
  inline virtual double F(double *);
  inline virtual double E(double *){assert(1); return 0;};
  virtual double M(double *){assert(1); return 0;};
  virtual void derivative(double *, double *){}
};

#endif



