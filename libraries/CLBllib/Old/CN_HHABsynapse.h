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



#ifndef CN_HHABSYNAPSE_H
#define CN_HHABSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define HHABPNO 7

#define HHABIVARNO 0

class HHABsynapse: public synapse
{

 public:

	 HHABsynapse(neuron *, neuron *, double, double, double, double, double,
			double, double, int, int, int);
	HHABsynapse(neuron *, neuron *, double, double, double, double, double,
			double, double);
	HHABsynapse(neuron *, neuron *, double *);

	virtual ~HHABsynapse();
	virtual double gsyn();
	virtual void set_gsyn(double);
	virtual double Isyn(double *);
	virtual void derivative(double *, double *) {
	}
};

#endif
