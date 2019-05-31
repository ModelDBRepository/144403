/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_S01SYNAPSEECPLAST3_CC
#define CN_S01SYNAPSEECPLAST3_CC

#include "CN_S01synapse.cc"

#define TINFIN 1e20

S01synapseECplast3::S01synapseECplast3(neuron *insource, neuron *intarget,
				   double inksyn, double inEsyn, double inEpre,
				   double inasyn, double inbsyn,
				   double inVslope, 
				   double inlrnampl, double indelayT):
  S01synapse(insource, intarget, inksyn, inEsyn, inEpre, inasyn,
	    inbsyn, inVslope,
	    S01ECPLAST3IVARNO, S01ECPLAST3PNO, S01ECPLAST3)
{
  p[6]= inlrnampl;
  p[7]= indelayT;
  nextT= TINFIN;
  synapse_change= 0;
} 

S01synapseECplast3::~S01synapseECplast3()
{
}

void S01synapseECplast3::update_gsyn(double *x)
{
  static double chng; 
  if (x[0] > nextT) {  // x[0] is the time ...
    p[0]+= nextdg;
    if (p[0] < 0.0) p[0]= 0.0;
    if (!chngTq.empty()) {
      nextT= chngTq.pop();
      nextdg= dgq.pop();
    }
    else {
      nextT= TINFIN;
    }
    synapse_change= 1;
  }
  else synapse_change= 0;

  if ((source->start_spiking) || (target->start_spiking)) {
    if ((source->spike_time > 0.0) && (target->spike_time > 0.0)) {
      double tau= target->spike_time-source->spike_time;
      chng= STDP_func(tau);
      if (chng != 0.0) {
	if (nextT == TINFIN) {
	  nextT= x[0]+p[7];
	  nextdg= p[6]*chng;
	}
	else {
	  chngTq.push(x[0]+p[7]);
	  dgq.push(p[6]*chng);
	}
      }
    }
  }
}

#define SHIFT 0.0
#define a1 -0.000000000413071
#define a2 1.632662605846178
#define a3 0.000000000031599
#define a4 -1.398394710295993

double S01synapseECplast3::STDP_func(double t)
{
  if (t+SHIFT < 0.0) {
    return a1*pow(t+SHIFT,16)*exp((t+SHIFT)*a2);
  }
  else {
    return a3*pow(t+SHIFT,16)*exp((t+SHIFT)*a4);
  }
}

#undef SHIFT
#undef a1
#undef a2
#undef a3
#undef a4

// end of class implementation

#undef TINFIN

#endif


