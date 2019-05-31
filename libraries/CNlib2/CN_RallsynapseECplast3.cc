/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_RALLSYNAPSEECPLAST3_CC
#define CN_RALLSYNAPSEECPLAST3_CC

#include "CN_Rallsynapse.cc"

#define TINFIN 1e20

RallsynapseECplast3::RallsynapseECplast3(neuron *insource, neuron *intarget,
				   double ingsyn, double inEsyn, double inEpre,
				   double intau, double inlrnampl,
					 double indelayT): 
  Rallsynapse(insource, intarget, ingsyn, inEsyn, inEpre, intau,
	    RALLECPLAST3IVARNO, RALLECPLAST3PNO, RALLECPLAST3)
{
  p[4]= inlrnampl;
  p[5]= indelayT;
  nextT= TINFIN;
  synapse_change= 0;
} 

RallsynapseECplast3::~RallsynapseECplast3()
{
}

void RallsynapseECplast3::update_gsyn(double *x)
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
	  nextT= x[0]+p[5];
	  nextdg= p[4]*chng;
	}
	else {
	  chngTq.push(x[0]+p[5]);
	  dgq.push(p[4]*chng);
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

#define SHIFT 0.0
#define a1 -0.000000000413071
#define a2 1.632662605846178
#define a3 0.000000000031599
#define a4 -1.398394710295993

double RallsynapseECplast3::STDP_func(double t)
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


