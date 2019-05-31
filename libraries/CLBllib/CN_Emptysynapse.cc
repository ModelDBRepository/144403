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

#ifndef CN_EMPTYSYNAPSE_CC
#define CN_EMPTYSYNAPSE_CC

#include "CN_synapse.cc"

// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag

Emptysynapse::Emptysynapse(neuron *insource, neuron *intarget, double ingsyn,
		       int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
}

// This is the constructor to be used directly ...

Emptysynapse::Emptysynapse(neuron *insource, neuron *intarget, double ingsyn):
  synapse(insource, intarget, EMPTYSYNAPSEIVARNO, EMPTYSYNAPSEPNO, EMPTYSYNAPSE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
}

Emptysynapse::Emptysynapse(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, EMPTYSYNAPSEIVARNO, EMPTYSYNAPSEPNO, EMPTYSYNAPSE)
{
  set_p(inp);
}

Emptysynapse::~Emptysynapse()
{
}

double Emptysynapse::gsyn()
{
  return p[0];
}

void Emptysynapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double Emptysynapse::Isyn(double *x)
{
//	return p[0]*source->S(x)*(vrest-source->p[12]);
  return p[0]*source->S(x)*(target->E(x)-source->p[12]);
}

// end of class implementation

#endif

