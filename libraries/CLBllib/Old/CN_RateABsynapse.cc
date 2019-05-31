   /*--------------------------------------------------------------------------
   Author: Thomas Nowotny

   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402

   email to:  tnowotny@ucsd.edu

   initial version: 2005-08-17

--------------------------------------------------------------------------*/

#ifndef CN_RATEABSYNAPSE_CC
#define CN_RATEABSYNAPSE_CC

#include "CN_synapse.cc"


// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag

RateABsynapse::RateABsynapse(neuron *insource, neuron *intarget,
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

RateABsynapse::RateABsynapse(neuron *insource, neuron *intarget,
		     double ingsyn, double inEsyn, double inEpre,
		     double inasyn, double inbsyn, double inrtime,
		     double nslevel):
  synapse(insource, intarget, RATEABSYNIVARNO, RATEABSYNPNO, RATEABSYN)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= inasyn;           // alpha timescale in 1/msec
  p[4]= inbsyn;           // beta timescale in 1/msec
  p[5]= inrtime;       // time of transmitter release
  p[6]= nslevel;          // noise level
}

RateABsynapse::RateABsynapse(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, RATEABSYNIVARNO, RATEABSYNPNO, RATEABSYN)
{
  set_p(inp);
}

RateABsynapse::~RateABsynapse()
{
}

double RateABsynapse::gsyn()
{
  return p[0];
}

void RateABsynapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double RateABsynapse::S(double *x)
{
  assert(enabled);
  return x[idx];
}


double RateABsynapse::Isyn(double *x)
{
  double vrest = -63.4675;
  double vrev = p[1];
  return -p[0]*x[idx]*(vrest-vrev);
}


void RateABsynapse::derivative(double *x, double *dx)
{
  double alpha = p[3];
  double beta = p[4];
  double tr = p[5];
  dx[idx] =  -x[idx]*beta+(alpha)*(exp(beta*tr)-1.0)/(exp(beta/source->F(x))-1);

}

void RateABsynapse::noise(double *x, double *dx)
{
 // dx[idx]= p[6]*abs(RG.n());
}

void RateABsynapse::init(double *x, double *iniVars)
{
  synapse::init(x, iniVars);
}

// end of class implementation

#endif


