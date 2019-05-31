/*--------------------------------------------------------------------------
   Author: Christopher Buckley



--------------------------------------------------------------------------*/

#ifndef CN_POISSONRATENEURON_CC
#define CN_POISSONRATENEURON_CC

#include "CN_neuron.cc"

PoissonRateNeuron::PoissonRateNeuron(int inlabel, double *the_p= POISSONRATENEURON_p):
  neuron(inlabel, POISSONRATENEURON_IVARNO, POISSONRATENEURON, the_p, NEUROSYN_PNO)
{
	mVrest = -60;
}

PoissonRateNeuron::PoissonRateNeuron(int inlabel, vector<int> inpos, double *the_p= POISSONRATENEURON_p):
  neuron(inlabel, POISSONRATENEURON_IVARNO, POISSONRATENEURON, inpos, the_p, NEUROSYN_PNO)
{
}

inline double PoissonRateNeuron::E(double *x)
{
  assert(enabled);
  return mVrest;
}


inline double PoissonRateNeuron::F(double *x)
{
	return Isyn;
}

double PoissonRateNeuron::S(double *x)
{
  assert(enabled);
  return x[idx];
}


void PoissonRateNeuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }

  double alpha = p[5];
  double beta = p[6];
  double tr = p[7];

  dx[idx] =  -x[idx]*beta+alpha*tr*Isyn;
}

void PoissonRateNeuron::init(double *x, double *iniVars)
{
	assert(enabled);
	  for (int i= 0; i < iVarNo; i++)
	    x[idx+i]= 0.0;
}



void PoissonRateNeuron::ResetSynapse(double *x){
	x[idx] = 0;
}

#endif
