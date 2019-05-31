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

#ifndef CN_HHABSYNAPSE_CC
#define CN_HHABSYNAPSE_CC

#include "CN_synapse.cc"

// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag


HHABsynapse::HHABsynapse(neuron *insource, neuron *intarget,
		     double ingsyn, double inEsyn, double inEpre,
		     double inasyn, double inbsyn, double inrtime,
		     double nslevel,
		     int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= inasyn;           // alpha timescale in 1/msec
  p[4]= inbsyn;           // beta timescale in 1/msec
  p[5]= inrtime;       // time of transmitter release
  p[6]= nslevel;          // noise level
}



// This is the constructor to be used directly ...

HHABsynapse::HHABsynapse(neuron *insource, neuron *intarget,
		     double ingsyn, double inEsyn, double inEpre,
		     double inasyn, double inbsyn, double inrtime,
		     double nslevel):
  synapse(insource, intarget, HHABIVARNO, HHABPNO, HHABSYN)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= inasyn;           // alpha timescale in 1/msec
  p[4]= inbsyn;           // beta timescale in 1/msec
  p[5]= inrtime;       // time of transmitter release
  p[6]= nslevel;          // noise level
}



HHABsynapse::HHABsynapse(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, HHABIVARNO, HHABPNO, HHABSYN)
{
  set_p(inp);
}

HHABsynapse::~HHABsynapse()
{
}

double HHABsynapse::gsyn()
{
  return p[0];
}

void HHABsynapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double HHABsynapse::Isyn(double *x)
{

double vrest = -63.4675;
  return p[0]*(vrest-p[1])*source->E(x);
}

// end of class implementation

#endif


