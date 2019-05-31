/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_POPPOISSONN_CC
#define CN_POPPOISSONN_CC

#include "CN_neuron.cc"

PopPoissonN::PopPoissonN(int inlabel, vector<int> inpos,
			     double *the_p= POPPOI_p):
  neuron(inlabel, POPPOI_IVARNO, POPPOISSONN,
	      inpos, the_p, POPPOI_PNO)
{
  fire_t= -1e10;
  next_spike= -1;
  setIdx(-1);
}

PopPoissonN::PopPoissonN(int inlabel, double *the_p= POPPOI_p):
  neuron(inlabel, POPPOI_IVARNO, POPPOISSONN,
	 the_p, POPPOI_PNO)
{
  fire_t= -1e10;
  next_spike= -1;
  setIdx(-1);
}

PopPoissonN::~PopPoissonN()
{
}

double PopPoissonN::E(double *x)
{
  if (x[0] - fire_t < p[3]) {  // spiking
    return p[2];
  }
  else {
    return p[1];  // baseline value
  }
}

void PopPoissonN::validate_E(double *x, double dt)
{
  if (x[0] - fire_t < p[0]) {
    // do nothing we are still spiking/refractory
  }
  else {
    if (R.n() <= *Lambda*dt) {
      next_spike= x[0];
    }
  }
}

void PopPoissonN::step()
{
  if (next_spike > 0) {
    fire_t= next_spike;
    next_spike= -1;
  }
}

void PopPoissonN::init(double *x, double *iniVars)
{
  fire_t= -1e10;
  next_spike= -1e10;
}

// allow to reference an external parameter Lambda
// so that it does not have to be reset all the time
void PopPoissonN::set_Lambda(double *Lbd)
{
  Lambda= Lbd;
}

#endif





