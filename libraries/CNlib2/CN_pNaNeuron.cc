/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_PNANEURON_CC
#define CN_PNANEURON_CC

#include "CN_neuron.cc"

pNaNeuron::pNaNeuron(int inlabel, double *the_p= pNa_p):
  neuron(inlabel, pNa_IVARNO, PNANEURON, the_p, pNa_PNO)
{
}

pNaNeuron::pNaNeuron(int inlabel, vector<int> inpos, double *the_p= pNa_p):
  neuron(inlabel, pNa_IVARNO, PNANEURON, inpos, the_p, pNa_PNO)
{
}

inline double pNaNeuron::E(double *x)
{
  assert(enabled);
  return x[idx];
}

#define gpNa p[0]
#define ENa p[1]
#define VmpNa p[2]
#define smpNa p[3]
#define taumpNa p[4]
#define VhpNa p[5]
#define shpNa p[6]
#define ChpNa p[7]
#define tauh0pNa p[8]
#define tauhApNa p[9]
#define VthpNa p[10]
#define sthpNa p[11]
#define Cmem p[12]

#define _efunc(a,b,V) 1.0/(1.0 + exp(((V)-(a))/(b)))

void pNaNeuron::currents(ostream &os, double *x)
{
  os << -x[idx+1]*x[idx+2]*gpNa*(x[idx]-ENa);
}

    
void pNaNeuron::derivative(double *x, double *dx)
{
  static double minf, hinf, taum, tauh;

  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }
  
  // differential eqn for E, the membrane potential
  dx[idx]= -(x[idx+1]*x[idx+2]*gpNa*(x[idx]-ENa)-Isyn)/p[9];

  // differential equation for pNa current activation 
  minf= _efunc(VmpNa, smpNa, x[idx]);
  taum= taumpNa; // constant time scale for activation (?)
  dx[idx+1]= (minf - x[idx+1])/taum;

  // differential equation for pNa current inactivation 
  hinf= (1.0-ChpNa)*(_efunc(VhpNa, shpNa, x[idx]))+ChpNa;
  tauh= tauh0pNa - tauhApNa*_efunc(VthpNa, sthpNa, x[idx]); // wild guess
  dx[idx+2]= (hinf - x[idx+2])/tauh;
}

#undef _efunc
#undef gpNa
#undef ENa
#undef VmpNa
#undef smpNa
#undef taumpNa
#undef VhpNa
#undef shpNa
#undef ChpNa
#undef tauh0pNa
#undef tauhApNa
#undef VthpNa
#undef sthpNa
#undef Cmem

#endif
