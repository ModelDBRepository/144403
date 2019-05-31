/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSE_SMSTDP1_CC
#define CN_ABSYNAPSE_SMSTDP1_CC

#include "CN_absynapse_smSTDP.cc"

// This is the constructor to be used directly ...

absynapse_smSTDP1::absynapse_smSTDP1(neuron *insource, neuron *intarget,
				     double inEsyn, double inEpre,
				     double inasyn,
				     double inbsyn, double inVslope,
				     double indecay, double ingmax,
				     double ing0, double ingdecay,
				     double inAp, double inAm,
				     double intaup, double intaum
				     ):
  absynapse_smSTDP(insource, intarget, SMSTDPIVARNO, SMSTDP1PNO, SMSTDP1)
{
  p[0]= inEsyn;           // Esyn reversal potential in mV
  p[1]= inEpre;           // Epre presyn threshold potential in mV
  p[2]= inasyn;           // alpha timescale in 1/msec
  p[3]= inbsyn;           // beta timescale in 1/msec
  p[4]= inVslope;         // steepness of activation curve as func of Vpre
  p[5]= indecay;
  p[6]= ingmax;
  p[7]= ing0;
  p[8]= ingdecay;
  p[9]= inAp;
  p[10]= inAm;
  p[11]= intaup;
  p[12]= intaum;
} 

absynapse_smSTDP1::absynapse_smSTDP1(neuron *insource, neuron *intarget, double *inp):
  absynapse_smSTDP(insource, intarget, SMSTDPIVARNO, SMSTDP1PNO, SMSTDP1)
{
  set_p(inp);
} 

absynapse_smSTDP1::~absynapse_smSTDP1()
{
}

// simple fit to Bi and Poo

double absynapse_smSTDP1::stdp_fn(double dt)
{
  if (dt > 0.0) {
    return p[9]*dt*exp(-dt/p[11]);
  }
  else {
    return p[10]*dt*exp(dt/p[12]);
  }
}

#endif


