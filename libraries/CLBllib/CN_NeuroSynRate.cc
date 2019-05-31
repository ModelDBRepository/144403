/*--------------------------------------------------------------------------
   Author: Thomas Nowotny

   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402

   email to:  tnowotny@ucsd.edu

   initial version: 2005-08-17

--------------------------------------------------------------------------*/

#ifndef CN_NEUROSYNRATE_CC
#define CN_NEUROSYNRATE_CC

#include "CN_neuron.cc"

NeuroSynRate::NeuroSynRate(int inlabel, double *the_p= NEUROSYNRATE_p):
  neuron(inlabel, NEUROSYNRATE_IVARNO, NEUROSYNRATE, the_p, NEUROSYN_PNO)
{
	mVrest = -55;;
}

NeuroSynRate::NeuroSynRate(int inlabel, vector<int> inpos, double *the_p= NEUROSYNRATE_p):
  neuron(inlabel, NEUROSYNRATE_IVARNO, NEUROSYNRATE, inpos, the_p, NEUROSYN_PNO)
{
}

inline double NeuroSynRate::E(double *x)
{
  assert(enabled);
  return mVrest;
}


inline double NeuroSynRate::F(double *x)
{
	return mF;
}

double NeuroSynRate::S(double *x)
{
  assert(enabled);
  return x[idx];
}


void NeuroSynRate::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }

  double beta = p[15];

  dx[idx] =  -x[idx]*beta+Transfer(Isyn+p[11]);
}

void NeuroSynRate::init(double *x, double *iniVars)
{
	assert(enabled);
	mF=0;
	  for (int i= 0; i < iVarNo; i++)
	    x[idx+i]= iniVars[i];
}



double NeuroSynRate::Transfer(double isyn) {

	double alpha = p[14];
	double tr = p[16];
	double currentIn;

	currentIn= isyn;
	//currentIn = max(isyn, 0.0);


		mF = M * currentIn + C;
		mF = max(mF, 0.0);
		mF = min(mF, 0.2);

	//	if(mF<0.004)
		//	mF = 0.0;
		return alpha * tr * mF;
}



void NeuroSynRate::ResetSynapse(double *x){
	x[idx] = 0;
}

#endif
