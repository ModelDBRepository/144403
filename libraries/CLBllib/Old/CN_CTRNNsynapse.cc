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

#ifndef CN_CTRNNSYNAPSE_CC
#define CN_CTRNNSYNAPSE_CC

#include "CN_synapse.cc"

// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag

CTRNNsynapse::CTRNNsynapse(neuron *insource, neuron *intarget, double ingsyn,
		       int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
}

// This is the constructor to be used directly ...

CTRNNsynapse::CTRNNsynapse(neuron *insource, neuron *intarget, double ingsyn):
  synapse(insource, intarget, CTRNNIVARNO, CTRNNPNO, CTRNNSYN)
{
  p[0]= ingsyn;           // gsyn strength of synapse
}

CTRNNsynapse::CTRNNsynapse(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, CTRNNIVARNO, CTRNNPNO, CTRNNSYN)
{
  set_p(inp);
}

CTRNNsynapse::~CTRNNsynapse()
{
}

double CTRNNsynapse::gsyn()
{
  return p[0];
}

void CTRNNsynapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double CTRNNsynapse::Isyn(double *x)
{
  return p[0]*source->E(x);
}

// end of class implementation

#endif


