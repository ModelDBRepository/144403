/*--------------------------------------------------------------------------
   Author: Thomas Nowotny

   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402

   email to:  tnowotny@ucsd.edu

   initial version: 2005-08-17

--------------------------------------------------------------------------*/

#ifndef CN_RATEPNNEURON_CC
#define CN_RATEPNNEURON_CC

#include "CN_neuron.cc"

RatePNneuron::RatePNneuron(int inlabel, double *the_p= RATEPN_p):
  neuron(inlabel, RATEPN_IVARNO, RATEPNNEURON, the_p, RATEPN_PNO)
{
	enabled =1;
}

RatePNneuron::RatePNneuron(int inlabel, vector<int> inpos, double *the_p= RATEPN_p):
  neuron(inlabel, RATEPN_IVARNO, RATEPNNEURON, inpos, the_p, RATEPN_PNO)
{
	enabled =1;
}

inline double RatePNneuron::F(double *x)
{

	Isyn= 0.0;
	forall(den, den_it) {
	    Isyn+= (*den_it)->Isyn(x);
	  }

  double A= 0.186337119371455;
  double betaF =  0.7176;
  double vrest = -63.4675;
  double Ic= -p[11]+0.041270138525220;;
  double g = p[10];
  double Ek = p[3];
  double D= vrest - Ek;

  double F;


  double A2 = pow(A,2);

  double currentIn = max(Isyn-Ic,0.0);
  F= (sqrt(4*A2*(currentIn) + pow((A2*betaF*g*D),2)) - A2*betaF*g*D)/2;


  return max(F,0.0);





}





#endif
