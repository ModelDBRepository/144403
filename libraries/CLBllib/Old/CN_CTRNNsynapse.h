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



#ifndef CN_CTRNNSYNAPSE_H
#define CN_CTRNNSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define CTRNNPNO 1

#define CTRNNIVARNO 0

class CTRNNsynapse: public synapse
{

 public:
  CTRNNsynapse(neuron *, neuron *, double, int, int, int);
  CTRNNsynapse(neuron *, neuron *, double);
  CTRNNsynapse(neuron *, neuron *, double *);
  virtual ~CTRNNsynapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *) { }
};

#endif
