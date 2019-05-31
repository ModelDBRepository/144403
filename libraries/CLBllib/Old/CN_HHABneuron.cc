//--------------------------------------------------------------------------
// Author: Christopher L. Buckley
//
// Institute: CCNR
//			  Informatics
//            University of Sussex
//
// email to:  c.l.buckley@sussex.ac.uk
//
// initial version: 2009-12-15
//
//--------------------------------------------------------------------------


#ifndef CN_HHABNEURON_CC
#define CN_HHABNEURON_CC


#include "CN_neuron.cc"

HHABneuron::HHABneuron(int inlabel, double *the_p= HHAB_p):
  neuron(inlabel, HHAB_IVARNO, HHABNEURON, the_p, HHAB_PNO)
{
}

HHABneuron::HHABneuron(int inlabel, vector<int> inpos, double *the_p= HHAB_p):
  neuron(inlabel, HHAB_IVARNO, HHABNEURON, inpos, the_p, HHAB_PNO)
{
}

void HHABneuron::set_p(double *inpn,double *inps)
{
  for (int i= 0; i < HHAB_NPNO; i++) p[i]= inpn[i];
  for (int i= 0; i < HHAB_SPNO; i++) {p[i+HHAB_NPNO]= inps[i];}
}

inline double HHABneuron::E(double *x)
{
  assert(enabled);
  return x[idx];
}

void HHABneuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }

  // differential eqn for E, the membrane potential





  double A= 0.186337119371455;
  double betaF =  0.7176;
  double Ic= p[11];
  double alpha = p[15];
  double beta = p[16];
  double tr = p[17];
  double g = p[10];
  double vrest = -63.4675;
  double vreverse = p[13];
  double D= vrest - vreverse;




  double F;


  double A2 = pow(A,2);

  F= (sqrt(4*A2*(Isyn-Ic) + pow((A2*betaF*g*D),2)) - A2*betaF*g*D)/2;

  dx[idx]= -beta*x[idx] + (alpha *exp(beta*tr)-1.0)/(exp(beta/F)-1);

}

void HHABneuron::noise(double *x, double *dx)
{
  // differential eqn for E, the membrane potential
 // dx[idx]= p[1]*RG.n();
}

#endif
