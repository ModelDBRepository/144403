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



#ifndef CN_EMPTYSYNAPSE_H
#define CN_EMPTYSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define EMPTYSYNAPSEPNO 1

#define EMPTYSYNAPSEIVARNO 0

class Emptysynapse: public synapse
{

 public:
  Emptysynapse(neuron *, neuron *, double, int, int, int);
  Emptysynapse(neuron *, neuron *, double);
  Emptysynapse(neuron *, neuron *, double *);
  virtual ~Emptysynapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *) { }
};

#endif
