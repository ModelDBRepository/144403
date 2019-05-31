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


#ifndef CN_PNANEURON_H
#define CN_PNANEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of persitent sodium current in Lymnaea cells

#define pNa_IVARNO 3
#define pNa_PNO 13

double stdpNa_p[pNa_PNO]= {
  0.233,          // 0 - gpNa: max conductance of persistent Na [1/(mOhms * cm^2)]
  30.0,          // 1 - ENa: Na equi potential [mV]
  -18,          // 2 - VmpNa: activation mid point [mV]
  -15.0,         // 3 - smpNa: slope of activation [mV]
  5.0,          // 4 - taumpNa: time scale of activation m [ms]
  -46,         // 5 - VhpNa: inactivation mid point [mV]
  7.43,         // 6 - shpNa: slope of inactivation [mV]
  0.11,         // 7 - ChpNa: persistent part of inactivation [unitless]
  50.0, //125        // 8 - tauh0pNa: maximum time scale of inactivation h [ms]
  25.0, //100        // 9 - tauhApNa: amplitude of tauh voltage dependence [ms]
  -70.0,         // 10 - VthpNa: mid point of tauh curve [mV]
  10.0,          // 11 - sthpNa: slope parameter of tauh curve [mV]
  1.0,         // 12 - Cmem: membr. capacity density in [muF/cm^2]
};

double *pNa_p= stdpNa_p;

const char *pNa_p_text[pNa_PNO]= {
  "0 - gpNa: max conductance of persistent Na [1/(mOhms * cm^2)]",
  "1 - ENa: Na equi potential [mV]",
  "2 - VmpNa: activation mid point [mV]",
  "3 - smpNa: slope of activation [mV]",
  "4 - taumpNa: time scale of activation m [ms]",
  "5 - VhpNa: inactivation mid point [mV]",
  "6 - shpNa: slope of inactivation [mV]",
  "7 - ChpNa: persistent part of inactivation [unitless]",
  "8 - tauh0pNa: maximum time scale of inactivation h [ms]",
  "9 - tauhApNa: amplitude of tauh voltage dependence [ms]",
  "10 - VthpNa: mid point of tauh curve [mV]",
  "11 - sthpNa: slope parameter of tauh curve [mV]",
  "12 - Cmem: membr. capacity density in [muF/cm^2]"
};

double pNa_INIVARS[pNa_IVARNO]= {
  -64.1251,                       // 0 - membrane potential E
  0.0176331,                   // 1 - prob. for Na channel activation m
  0.994931,                   // 2 - prob. for not Na channel blocking h
};

const char *pNa_INIVARSTEXT[pNa_IVARNO]= {
  "0 - membrane potential E",
  "1 - pNa channel activation mpNa",
  "2 - pNa channel inactivation hpNa",
};


// PNaNeuron class itself

class pNaNeuron: public neuron
{
 private:
  double Isyn;

 public:
  pNaNeuron(int, double *);
  pNaNeuron(int, vector<int>, double *);
  ~pNaNeuron() { }
  inline virtual double E(double *);
  virtual void currents(ostream &, double *);
  virtual void derivative(double *, double *);
};

#endif



