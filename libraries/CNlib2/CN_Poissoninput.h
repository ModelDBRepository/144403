//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institute for Nonlinear Dynamics
//            University of California San Diego
//            La Jolla, CA 92093-0402
//
// email to:  tnowotny@ucsd.edu
//
// initial version: 2003-06-18
//
//--------------------------------------------------------------------------


#ifndef POISSONINPUT_H
#define POISSONINPUT_H

#include <cmath>
#include "randomGen.h"
#include "randomGen.cc"

randomGen R;
#define POISSONINPUT 101

#define POI_IVARNO 0
#define POI_PNO 18

double stdPOI_p[POI_PNO]= {
  0.1,           // 0 - Lambda: firing rate
  0.0,           // 1 - refratory period
  20.0,          // 2 - Vspike
  -60.0 ,         // 3 - Vrest
  0.021, // 4 - gl: leak conductance in 1/(mOhms * cm^2)
  		-55.0, // 5 - El: leak equi potential in mV
  		0.00572, // 6 - gKl: potassium leakage conductivity
  		-95.0, // 7 - EKl: potassium leakage equi pot in mV
  		65.0, // 8 - V0: ~ total equi potential (?)
  		0.143, // 9 - Cmem: membr. capacity density in muF/cm^2
  		0,//0.715, // 10 - gM: conductance of the M current
  		0.0, // 11- IDC: baseline offset current
  		-60, //  12 inEsyn %reversal potential (-95 = inhibitory)
  		-20, // 13  inEpre %threshold for pre synaptic spike detection
  		2, //14 inasyn
  				0.05, // 15 inbsyn
  				1, // 16 inrtime
  		0 //17 noise
};

double *POI_p= stdPOI_p;

const char *POI_p_text[POI_PNO]= {
  "0 - Lmabda: firing rate",
  "1 - refractory period",
  "2 - Vspike",
  "3 - Vrest"
};

// the POI neuron class itself

class Poissoninput:public neuron
{
 public:
  int Evalid;
  int Isynvalid;
  int refract;
  double mS,mSLast;

  Poissoninput(int, double *);
  ~Poissoninput();
  void set_input(double);
  void init();
  void advance(double *, double);
  double E(double *);
  void derivative(double *, double *) { }

  inline virtual double F(double *){return 0;}
  double S(double *);
};

#endif
