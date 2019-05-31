/*--------------------------------------------------------------------------
 Author: Thomas Nowotny

 Institute: Institute for Nonlinear Dynamics
 University of California San Diego
 La Jolla, CA 92093-0402

 email to:  tnowotny@ucsd.edu

 initial version: 2005-08-17

 --------------------------------------------------------------------------*/

#ifndef CN_POISSONNEURON_H
#define CN_POISSONNEURON_H

#include "CN_neuron.h"

#define POI_IVARNO 1
double *POI_p = PARAMS_ORN;

class Poissonneuron: public neuron {
private:
	double fire_t;
	double tlast;
	int firing;
	int refract;
	double *myx, *myxn;
	double Isyn;

public:
	Poissonneuron(int, vector<int> , double *);
	Poissonneuron(int, double *);
	~Poissonneuron();
	 virtual double F(double *);
	double S(double *);
	virtual double E(double *);
	virtual void validate_E(double *, double);
	virtual void step();
	virtual void derivative(double *, double *);
	virtual void init(double *, double *);
	bool spiking(void);
	  double mVrest;
	  double spike;
};

#endif

