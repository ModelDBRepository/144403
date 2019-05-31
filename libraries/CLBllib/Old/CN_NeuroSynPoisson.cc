/*--------------------------------------------------------------------------
   Author: Thomas Nowotny

   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402

   email to:  tnowotny@ucsd.edu

   initial version: 2005-08-17

--------------------------------------------------------------------------*/

#ifndef CN_NEUROSYNPOISSON_CC
#define CN_NEUROSYNPOISSON_CC

#include "CN_neuron.cc"

NeuroSynPoisson::NeuroSynPoisson(int inlabel, double *the_p= NEUROSYNPOISSON_p):
  neuron(inlabel, NEUROSYNPOISSON_IVARNO, NEUROSYNPOISSON, the_p, NEUROSYNPOISSON_PNO)
{
	 tlast= -10000;
}

NeuroSynPoisson::NeuroSynPoisson(int inlabel, vector<int> inpos, double *the_p= NEUROSYNPOISSON_p):
  neuron(inlabel, NEUROSYNPOISSON_IVARNO, NEUROSYNPOISSON, inpos, the_p, NEUROSYNPOISSON_PNO)
{
	 tlast= -10000;
}



double NeuroSynPoisson::E(double *x)
{
  if (spikingThistime) {
    return p[20];
  }
  else {
    return p[21];
  }
}

double NeuroSynPoisson::S(double *x)
{
  assert(enabled);
  return x[idx];
}
void NeuroSynPoisson::derivative(double *x, double *dx)
{


  static double dt;

   dt= x[0] - tlast;
   if ((dt >= 0) && (dt <= p[16])) {
 	  dx[idx]= p[14] - p[15]*x[idx];
   } else {
     if ((E(x) > p[13]) && (dt > p[16])) {
       // new spike ... start releasing
       tlast= x[0];
      dx[idx]= p[14] - p[15]*x[idx];
     }
     else {
       // no release
       dx[idx]= -p[15]*x[idx];
     }
   }

   if (spiking)
       spiking= 0;

}

void NeuroSynPoisson::init(double *x, double *iniVars)
{

	assert(enabled);
	  for (int i= 0; i < iVarNo; i++)
	    x[idx+i]= iniVars[i];

	  spike_time= -1.0;
	    spiking= 0;
	    refract= 0;


	    spikingThistime = 0.0;

	  start_spiking= 0;
	  tlast= -10000;
}

void NeuroSynPoisson::ResetSynapse(double *x)
{
	x[idx] = 0;
}

void NeuroSynPoisson::advance(double *x, double dt,double buts)
{

	if (spikingThistime)
			   spikingThistime =0.0;

	 if (spiking)
	    spiking= 0;

	 double actualTime = x[0]-buts/p[18];

  if((actualTime - spike_time) > 1/p[18]){
	  spiking= 1;
	  spike_time= actualTime;
	  spikingThistime =1.0;
  }
/*

  if (RG.n() < p[18]*0.1) {
       spiking= 1;
       spikingThistime =1.0;
       spike_time= x[0];
     }
*/

}

#endif
