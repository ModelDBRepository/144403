/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-01-25
  
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSEECPLAST2_CC
#define CN_ABSYNAPSEECPLAST2_CC

#include "CN_absynapse.cc"

absynapseECplast2::absynapseECplast2(neuron *insource, neuron *intarget,
				   double inksyn, double inEsyn, double inEpre,
				   double inasyn, double inbsyn,
				   double inrtime, 
				   double inlrnampl):
  absynapse(insource, intarget, inksyn, inEsyn, inEpre, inasyn,
	    inbsyn, inrtime,
	    ABECPLAST2IVARNO, ABECPLAST2PNO, ABECPLAST2)
{
  p[6]= inlrnampl;
  synapse_change= 0;
} 

absynapseECplast2::~absynapseECplast2()
{
}

void absynapseECplast2::update_gsyn(double *x)
{
  static double chng; 
  if ((source->start_spiking) || (target->start_spiking)) {
    if ((source->spike_time > 0.0) && (target->spike_time > 0.0)) {
      double tau= target->spike_time-source->spike_time;
      chng= STDP_func(tau);
      if (chng != 0.0) {
	synapse_change= 1;
	p[0]+= p[6]*chng;
	if (p[0] < 0.0) p[0]= 0.0;
      }
      else synapse_change= 0;
    }
  }
}

#define SHIFT 5.0
#define a1 -0.000009954616079
#define a2 1.331251297482154
#define a3 0.000000212995136
#define a4 -0.896464339979975

double absynapseECplast2::STDP_func(double t)
{
//   if (t > 0.0) {  
//     return pow(t, 10.0)*exp(-abs(t))*5.14e-9;
//     // amplitude is fit to data ...
//   }
//   else {
//     return -pow(t, 10.0)*exp(-abs(t))*5.14e-9;
//   }
  if (t+SHIFT < 0.0) {
    return a1*pow(t+SHIFT,10)*exp((t+SHIFT)*a2);
  }
  else {
    return a3*pow(t+SHIFT,10)*exp((t+SHIFT)*a4);
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


