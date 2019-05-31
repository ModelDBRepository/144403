/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_DATA_CC
#define CN_DATA_CC

#include "CN_neuron.cc"

Data::Data(int inlabel, char *name, double *the_p= Data_p):
  neuron(inlabel, Data_IVARNO, DATA, the_p, Data_PNO)
{
  double data;
  ifstream *is= new ifstream(name);

  cerr << name << endl;
  assert(is->good());
  
  cnt= 0;
  while (!is->eof()) {
    *is >> data;
    cnt++;
  }
  dA= new double[cnt];
  delete is;

  is= new ifstream(name);
  for (int i= 0; i < cnt; i++) {
    *is >> dA[i];
  }
}

Data::~Data()
{
  delete[] dA;
}

inline double Data::E(double *x)
{
  static int pos;
  
  pos= (int) ((x[0]-p[0])/p[1]);
  if (pos < cnt) return dA[pos];
  else return 0.0;
}

#endif
