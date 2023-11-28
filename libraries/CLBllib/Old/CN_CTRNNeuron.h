//--------------------------------------------------------------------------
// Author: Christopher L. Buckley
//
// Institute: CCNR
//			  Informatics
//            University of Sussex
//
// email to:  c.l.buckley@sussex.ac.uk
//
// initial version: 2009-12-15
//
//--------------------------------------------------------------------------


#ifndef CN_CTRNNEURON_H
#define CN_CTRNNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define CTRNN_IVARNO 1
#define CTRNN_PNO 3

double stdCTRNN_p[CTRNN_PNO]= {
  1.0,          // 0 - rho_ii: "self inhibition"
  0.0,			 // bais value
  0.0           // 1 - sigma: noise level
};

double *CTRNN_p= stdCTRNN_p;

const char *CTRNN_p_text[CTRNN_PNO]= {
  "0 - rho_ii: self inhibition",
};

double CTRNN_INIVARS[CTRNN_IVARNO]= {
  0.1                       // 0 - firing rate
};

const char *CTRNN_INIVARSTEXT[CTRNN_IVARNO]= {
  "0 - firing rate"
};


// CTRNN neuron class itself

class CTRNNeuron: public neuron
{
 private:
  double Isyn;
 public:
  CTRNNeuron(int, double *);
  CTRNNeuron(int, vector<int>, double *);
  ~CTRNNeuron() { }
  inline virtual double E(double *);
  virtual void derivative(double *, double *);
  virtual void noise(double *, double *);
};

#endif



