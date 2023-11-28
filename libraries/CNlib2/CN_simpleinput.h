/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-02-14
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
  simple input (constant voltage)
--------------------------------------------------------------------------*/


#ifndef CN_SIMPLEINPUT_H
#define CN_SIMPLEINPUT_H

#include "CN_neuron.h"

#define SIN_IVARNO 1
#define SIN_PNO 1

double stdSIN_p[SIN_PNO]= {
  -60.0,          // 0 - standard potential
};

double *SIN_p= stdSIN_p;

const char *SIN_p_text[SIN_PNO]= {
  "0 - standard potential",
};

double SIN_INIVARS[SIN_IVARNO]= {
  -60.0           // 0 - potential
};

const char *SIN_INIVARSTEXT[SIN_IVARNO]= {
  "0 - potential"
};

class simpleinput: public neuron
{
 public:
  simpleinput(int, double *);
  simpleinput(int, vector<int>, double *);
  ~simpleinput();
  virtual void set_V(double *, double);
  virtual double E(double *);
  virtual void derivative(double *, double *) { }
};


#endif





