/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSEECPLAST1_CC
#define CN_ABSYNAPSEECPLAST1_CC

#include "CN_absynapse.cc"

absynapseECplast1::absynapseECplast1(neuron *insource, neuron *intarget,
				   double inksyn, double inEsyn, double inEpre,
				   double inasyn, double inbsyn,
				   double inrtime, double ingmax,
				   double ing12, double ingslope,
				   double inlrnampl):
  absynapse(insource, intarget, inksyn, inEsyn, inEpre, inasyn,
	    inbsyn, inrtime,
	    ABECPLAST1IVARNO, ABECPLAST1PNO, ABECPLAST1)
{
  p[6]= ingmax;           // g_max the maximum synapse strength
  p[7]= ing12;            // g_1/2 the midpoint of the smoothing fn
  p[8]= ingslope;         // g_slope the slope of the tanh
  p[9]= inlrnampl;
  set_gsyn(inksyn);       // ksyn strength of synapse
  synapse_change= 0;
} 

absynapseECplast1::~absynapseECplast1()
{
}

double absynapseECplast::rgsyn()
{
  // return raw gsyn
  return p[0];
}

void absynapseECplast::set_rgsyn(double ingsyn)
{
  // set raw gsyn
  p[0]= ingsyn;
}

double absynapseECplast::gsyn()
{
  double tmp= p[6]/2.0*(tanh(p[8]*(p[0]-p[7]))+1.0);
  return tmp;
}

void absynapseECplast::set_gsyn(double ingsyn)
{
  double tmp= ingsyn/p[6]*2.0-1.0;
  p[0]= 0.5*log((1.0+tmp)/(1.0-tmp))/p[8]+p[7];
}

double absynapseECplast::Isyn(double *x)
{
  return -gsyn()*iVars[0]*(target->E(x)-p[1]);
}

void absynapseECplast1::update_gsyn(double *x)
{
  static double chng; 
  if ((source->start_spiking) || (target->start_spiking)) {
    if ((source->spike_time > 0.0) && (target->spike_time > 0.0)) {
      double tau= target->spike_time-source->spike_time;
      chng= STDP_func(tau);
      if (chng != 0.0) {
	synapse_change= 1;
	p[0]+= p[9]*chng;
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

double absynapseECplast1::STDP_func(double t)
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


