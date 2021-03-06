/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

using namespace std;

#include <iostream>

#include "gauss.h"
#include "randomGen.h"
#include "randomGen.cc"

randomGen RG(1,2,3);
randomGen R(1,2,3);

#include "CN_ValAdaptneuron.h"
#include "CN_absynapse.h"
#include "CN_NeuronModel.h"
#include "CN_rk65n.h"
#include "CN_DCInput.h"

#include "CN_ValAdaptneuron.cc"
#include "CN_absynapse.cc"
#include "CN_NeuronModel.cc"
#include "CN_DCInput.cc"

int main(int argc, char *argv[])
{
  list<neuron *> neurs;
  list<synapse *> syns;
  neuron *n;
  // synapse *s;
  synapse *ins;
  double *x, *xn, *tmp;
  int N;

  cerr << RG.n() << endl;
  double IDC= atof(argv[1]);
  vector<int> pos(3, 0);
  n= new ValAdaptneuron(1,pos);
  neurs.push_back(n);
  //  s= new absynapse(n, n2, 0.5, 0.0, -20.0, 0.1, 0.2, 5.0, 0.0);
  //  syns.append(s);
  ins= new DCInput(n, IDC);
  syns.push_back(ins);
  
  NeuronModel model(&neurs, &syns, N, cerr);
  x= new double[N];
  xn= new double[N];
  //n->init(x, ECN_INIVARS);
  n->init(x, ValA_INIVARS);
  //  s->init(x, ABSYN_INIVARS);
  rk65n machine(N, 0.1, 1e-8, 1e-12, 1e-6);

  double dt= 0.1;
  double dtx= 0.1;
  for (int i= 0; i < 100000; i++) {
    dtx= machine.integrate(x, xn, &model, dt);
    //    cerr << dtx << endl;
    dt= dtx;
    tmp= x; x= xn; xn= tmp;
    cout << x[0] << " " << n->E(x) << endl;
    if (n->start_spiking) cerr << x[0] << endl;
  }
  return 0;
}
    
