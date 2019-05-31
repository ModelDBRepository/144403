/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#include "CN_InputFunction.h"
#include "CN_synapse.cc"

InputFunction:: InputFunction(neuron *target, double inA, double inOm):
  synapse((neuron *) NULL, target, DCIVARNO, DCPNO, INPUTFUNCTION)
{
  p[0]= inA;               // injected current
  p[1]= inOm;
}


InputFunction::~InputFunction()
{
}

void InputFunction::set_A(double A)
{
  p[0]= A;
}

void InputFunction::set_Om(double Om)
{
  p[1]= Om*2.0*3.1415926536;
}


double InputFunction::Isyn(double *x)
{
  return -p[0]*(sin(p[1]*x[0])+1.0)*(1.0-target->E(x));
}


