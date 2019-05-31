/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSEECPLAST3_CC
#define CN_ABSYNAPSEECPLAST3_CC

#include "CN_absynapse.cc"

#define TINFIN 1e20

absynapseECplast3::absynapseECplast3(neuron *insource, neuron *intarget,
				     double inksyn, double inEsyn, double inEpre,
				     double inasyn, double inbsyn,
				     double inrtime, double nslevel,
				     double inlrnampl, double indelayT):
  absynapse(insource, intarget, inksyn, inEsyn, inEpre, inasyn,
	    inbsyn, inrtime, nslevel, ABECPLAST3IVARNO, ABECPLAST3PNO,
	    ABECPLAST3) 
{
  p[6]= inlrnampl;
  p[7]= indelayT;
  nextT= TINFIN;
  synapse_change= 0;
} 

absynapseECplast3::~absynapseECplast3()
{
}

void absynapseECplast3::update_gsyn(double *x)
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

// #define SHIFT 5.0
// #define a1 -0.000009954616079
// #define a2 1.331251297482154
// #define a3 0.000000212995136
// #define a4 -0.896464339979975

// #define SHIFT 6.0
// #define a1 -0.013397823405251
// #define a2 0.706022930241506
// #define a3 0.000950123022128
// #define a4 -0.370845166583127

#define a1 -0.000000260113763
#define a2 0.943381530515740
#define a3 0.000002289108353
#define a4 -1.099830823861115

double absynapseECplast3::STDP_func(double t)
{
//   if (t > 0.0) {  
//     return pow(t, 10.0)*exp(-abs(t))*5.14e-9;
//     // amplitude is fit to data ...
//   }
//   else {
//     return -pow(t, 10.0)*exp(-abs(t))*5.14e-9;
//   }
//   if (t+SHIFT < 0.0) {
//     return a1*pow(t+SHIFT,10)*exp((t+SHIFT)*a2);
//   }
//   else {
//     return a3*pow(t+SHIFT,10)*exp((t+SHIFT)*a4);
//   }
// }

  if (t < 0.0) {
    return a1*pow(t,10)*exp(t*a2);
  }
  else {
    return a3*pow(t,10)*exp(t*a4);
  }
}

#undef a1
#undef a2
#undef a3
#undef a4

// end of class implementation

#undef TINFIN

#endif


