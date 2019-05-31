/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#include "CN_InputFunctionNoise.h"
#include "CN_synapse.cc"

InputFunctionNoise:: InputFunctionNoise(neuron *target, double inA, double inOm,
			      double inNL):
  synapse((neuron *) NULL, target, DCIVARNO, DCPNO, INPUTFUNCTION)
{
  p[0]= inA;               // injected current
  p[1]= inOm;
  p[2]= inNL;
}


InputFunctionNoise::~InputFunctionNoise()
{
}

void InputFunctionNoise::set_A(double A)
{
  p[0]= A;
}

void InputFunctionNoise::set_Om(double Om)
{
  p[1]= Om*2.0*3.1415926536;
}

void InputFunctionNoise::set_Noise(double NL)
{
  p[2]= NL;
}

double InputFunctionNoise::Isyn(double *x)
{
  return -p[0]*(sin(p[1]*x[0])+1.0)*(1.0-target->E(x))+p[2]*RG.n();
}


