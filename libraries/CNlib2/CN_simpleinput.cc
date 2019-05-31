/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-02-14
  
--------------------------------------------------------------------------*/


#ifndef CN_SIMPLEINPUT_CC
#define CN_SIMPLEINPUT_CC

#include "CN_simpleinput.h"
#include "CN_neuron.cc"

simpleinput::simpleinput(int inlabel, double *the_p= SIN_p):
  neuron(inlabel, SIN_IVARNO, SIN, the_p, SIN_PNO)
{
}

simpleinput::simpleinput(int inlabel, vector<int> inpos,
			 double *the_p= SIN_p):
  neuron(inlabel, SIN_IVARNO, SIN, inpos, the_p, SIN_PNO)
{
}

simpleinput::~simpleinput()
{
}

void simpleinput::set_V(double *x, double V)
{
  x[idx]= V;
}

double simpleinput::E(double *x)
{
  return x[idx];
}


#endif
