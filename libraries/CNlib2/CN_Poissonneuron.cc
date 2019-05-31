/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_POISSONNEURON_CC
#define CN_POISSONNEURON_CC

#include "CN_neuron.cc"

Poissonneuron::Poissonneuron(int inlabel, vector<int> inpos,
			     double *the_p= POI_p):
  neuron(inlabel, POI_IVARNO, POISSONNEURON,
	      inpos, the_p, NEUROSYNADAPT_PNO)
{
  firing= 0;
  refract= 0;
  tlast= -10000;
  fire_t= 0;
  myx= new double[1];
  myxn= new double[1];
  myx[0]= p[2];
  myxn[0]= p[2];
  spike = 0;
//  setIdx(1);
}

Poissonneuron::Poissonneuron(int inlabel, double *the_p= POI_p):
  neuron(inlabel, POI_IVARNO, POISSONNEURON,
	 the_p, NEUROSYNADAPT_PNO)
{
  firing= 0;
  refract= 0;
  tlast= -10000;
  fire_t= 0;
  myx= new double[1];
  myxn= new double[1];
  myx[0]= p[2];
  myxn[0]= p[2];
  spike=0;
 // setIdx(-1);
}

Poissonneuron::~Poissonneuron()
{
  delete[] myx;
  delete[] myxn;
}


double Poissonneuron::S(double *x)
{
  assert(enabled);

  return x[idx];
}

double Poissonneuron::E(double *x)
{
  return myx[0];
}
double Poissonneuron::F(double *x)
{
  return Isyn;
}

void Poissonneuron::validate_E(double *x, double ddt)
{

	 Isyn= 0.0;
	  forall(den, den_it) {
	    Isyn+= (*den_it)->Isyn(x);
	  }


	  if (firing) {
	    if (x[0] - fire_t > p[0]) { // remember: x[0] is the time
	      firing= 0;
	      refract= 1;
	    }
	  }
	  else {
	    if (refract) {
	      if (x[0] - fire_t > p[1]) refract=0;
	    }
	    else {
	      if (RG.n() <= Isyn*ddt) {
		firing= 1;
		fire_t= x[0];
	      }
	    }
	  }
	  if (firing) myxn[0]= p[3];
	  else myxn[0]= p[2];
}

void Poissonneuron::step()
{
  myx[0]= myxn[0];
}

void Poissonneuron::init(double *x, double *iniVars)
{
  myx[0]= -60.0;
  myxn[0]=-60.0;
  x[idx] =0;
}
void Poissonneuron::derivative(double *x, double *dx)
{


	double alpha = p[5];
	  double beta = p[6];
	  double tr = p[7];

	  static double dt;

	   dt= x[0] - tlast;
	   if ((dt >= 0) && (dt <= p[7])) {
	 	  dx[idx]= alpha- beta*x[idx];
	   } else {
	     if ((myx[0] > p[8]) && (dt > p[7])) {
	       // new spike ... start releasing
	       tlast= x[0];
	      dx[idx]= alpha - beta*x[idx];
	     }
	     else {
	       // no release
	       dx[idx]= -beta*x[idx];
	     }
	   }

}


#endif





