/*--------------------------------------------------------------------------
   Author: Thomas Nowotny

   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402

   email to:  tnowotny@ucsd.edu

   initial version: 2002-01-25

--------------------------------------------------------------------------*/

#ifndef POISSONINPUT_CC
#define POISSONINPUT_CC

Poissoninput::Poissoninput(int inlabel, double *the_p= POI_p):
  neuron(inlabel, POI_IVARNO, POISSONINPUT, the_p, POI_PNO)
{
  spiking= 0;
  refract= 0;
  spike_time= -1.0;
  Evalid= 0;
  Isynvalid= 0;
  mS = 0;
  mSLast =0.0;
}

Poissoninput::~Poissoninput()
{
  delete[] p;
}


void Poissoninput::set_input(double inLambda)
{
  p[0]= inLambda;
}

void Poissoninput::init()
{
  spike_time= -1.0;
  spiking= 0;
  refract= 0;
}

void Poissoninput::advance(double *x, double dt)
{
  if (spiking) {
    spiking= 0;
    refract= 1;
  }

  if (refract) {
    if (x[0] - spike_time > p[1]) {
      refract= 0;
    }
  }

  if (!refract) {
    if (R.n() < p[0]*dt) {
      spiking= 1;
      spike_time= x[0];
    }
  }


  double tsls = x[0] - spike_time;

	if ((tsls >= 0) && (tsls <= p[16])) {
		mS = mSLast + dt * (p[14] - p[15] * mSLast);
	} else {
		mS = mSLast + dt * (-p[15] * mSLast);
	}

	mSLast = mS;
}

double Poissoninput::E(double *x)
{
  if (spiking) {
    return p[2];
  }
  else {
    return p[3];
  }
}



double Poissoninput::S(double *x)
{
   return mS;
}

#endif
