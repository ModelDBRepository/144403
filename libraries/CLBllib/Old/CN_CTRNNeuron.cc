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


#ifndef CN_CTRNNEURON_CC
#define CN_CTRNNEURON_CC


#include "CN_neuron.cc"

CTRNNeuron::CTRNNeuron(int inlabel, double *the_p= CTRNN_p):
  neuron(inlabel, CTRNN_IVARNO, CTRNNEURON, the_p, CTRNN_PNO)
{
}

CTRNNeuron::CTRNNeuron(int inlabel, vector<int> inpos, double *the_p= CTRNN_p):
  neuron(inlabel, CTRNN_IVARNO, CTRNNEURON, inpos, the_p, CTRNN_PNO)
{
}

inline double CTRNNeuron::E(double *x)
{
  assert(enabled);
  return x[idx];
}

void CTRNNeuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }

  // differential eqn for E, the membrane potential
  dx[idx]= -x[idx] +tanh(p[0]*x[idx]+p[1]+ Isyn);

}

void CTRNNeuron::noise(double *x, double *dx)
{
  // differential eqn for E, the membrane potential
  dx[idx]= p[1]*RG.n();
}

#endif
