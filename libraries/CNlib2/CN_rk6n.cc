/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_RK6N_CC
#define CN_RK6N_CC

#include "CN_rk6n.h"
#include <iostream>

rk6n::rk6n(int dim)
{
  N= dim;
  for (int i= 0; i < 9; i++) {
    Y[i]= new double[N];
    F[i]= new double[N];
    for (int j= i; j < 8; j++) {
      a[i][j]= 0.0;
    }
  }
  
  a[1][0]= 0.111111111111111;
  a[2][0]= 0.555555555555556e-1;
  a[2][1]= 0.555555555555556e-1;
  a[3][0]= 0.416666666666667e-1;
  a[3][1]= 0.0;
  a[3][2]= 0.125;
  a[4][0]= 0.166666666666667;
  a[4][1]= 0.0;
  a[4][2]= -0.5;
  a[4][3]= 0.666666666666667;
  a[5][0]= 0.1875e+1;
  a[5][1]= 0.0;
  a[5][2]= -0.7875e+1;
  a[5][3]= 0.7e+1;
  a[5][4]= -0.5;
  a[6][0]= -0.4227272727272727e+1;
  a[6][1]= 0.0;
  a[6][2]= 0.176995738636364e+2;
  a[6][3]= -0.142883522727273e+2;
  a[6][4]= 0.522017045454545;
  a[6][5]= 0.104403409090909e+1;
  a[7][0]= 0.840622673179752e+1;
  a[7][1]= 0.0;
  a[7][2]= -0.337303717185049e+2;
  a[7][3]= 0.271460231129622e+2;
  a[7][4]= 0.342046929709216;
  a[7][5]= -0.184653767923258e+1;
  a[7][6]= 0.577349465373733;
  a[8][0]= 0.128104575163399;
  a[8][1]= 0.0;
  a[8][2]= 0.0;
  a[8][3]= -0.108433734939759;
  a[8][4]= 0.669375;
  a[8][5]= -0.146666666666667;
  a[8][6]= 0.284444444444444;
  a[8][7]= 0.173176381998583;

  b[0]= 0.567119155354449e-1;
  b[1]= 0.0;
  b[2]= 0.0;
  b[3]= 0.210909572355356;
  b[4]= 0.141490384615385;
  b[5]= 0.202051282051282;
  b[6]= 0.253186813186813;
  b[7]= 0.843679809736684e-1;
  b[8]= 0.512820512820513e-1;
}

rk6n::~rk6n()
{
  for (int i= 0; i < 9; i++) {
    delete[] Y[i];
    delete[] F[i];
  }
}

void rk6n::integrate(double *y, double *yn,
			NeuronModel *model, double dt)
{
  // calculate iterative terms rk6_Y[__i] and rk6_F[__i] (to sixth order)
  //  yn[0]= y[0]+dt;
  for (k= 0; k < N; k++) {
    Y[0][k]= y[k];
  }
  model->derivative(Y[0], F[0]);
  for (i= 1; i < 9; i++)
  {
    for (k= 0; k < N; k++)
    {
      aF= 0.0;
      for (j= 0; j < i; j++)
	aF+= a[i][j]*F[j][k];
      Y[i][k]= y[k]+dt*aF;
    }
    model->derivative(Y[i], F[i]);
  }

  // no build of 5th order scheme

  // sum up rk6_Y[__i] and rk6_F[__i] to build 6th order scheme -> yn
  for (k= 0; k < N; k++)
  {
    aF= 0.0;
    for (j= 0; j < 9; j++) aF+= b[j]*F[j][k];
    yn[k]= y[k]+ dt*aF;
  }
}

#endif
