/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-01-25
  
--------------------------------------------------------------------------*/

#ifndef LPRNEURON_CC
#define LPRNEURON_CC

#include "CN_neuron.cc"

#define efunc(X,Y,Z) (1.0/(1.0+exp(((X)-(Y))/(Z))))
#define Eaxon x[idx]
#define Esoma x[idx+15]

LPRneuron::LPRneuron(int inlabel, double *the_p= LPR_p):
  neuron(inlabel, LPR_IVARNO, LPRNEURON, the_p, LPR_PNO)
{
}

LPRneuron::LPRneuron(int inlabel, vector<int> inpos, double *the_p= LPR_p):
  neuron(inlabel, LPR_IVARNO, LPRNEURON, inpos, the_p, LPR_PNO)
{
}

inline double LPRneuron::E(double *x)
{
  using namespace LPR5;
  return x[idx+15]+p[V_shift];
}

void LPRneuron::currents(ostream &os, double *x)
{
  using namespace LPR5;
  static double INa, ICa, IoCa, Id, IA, Ih, Il, IM, IVV, tmp;

  // differential eqn for the axon membrane potential
  INa= pw3(x[idx+1])*x[idx+2]*p[g_Na]*(Eaxon-p[V_Na]);
  Id= pw4(x[idx+8])*p[g_Kd]*(Esoma-p[V_K]);
  IM= p[g_M]*x[idx+14]*(Eaxon-p[V_M]);
  Il= p[g_leak]*(Eaxon-p[V_leak]);
  IVV= p[g_VV]*(Esoma-Eaxon);
  
  os << INa << " ";
  os << Id << " ";
  os << IM << " ";
  os << Il << " ";
  os << IVV << " ";

  tmp= exp(Esoma/p[RTF]);
  ICa= (x[idx+3]*x[idx+4]*p[g_CaT]+x[idx+5]*p[g_CaS])*
    p[P_Ca]*Esoma*(x[idx+13]*tmp-p[Ca_out])/(tmp-1.0);
  IoCa= x[idx+6]*x[idx+7]*p[g_KCa]*(Esoma-p[V_K]);
  tmp= efunc(Esoma,p[V_aA],p[s_aA]);
  IA= pw3(x[idx+9])*p[g_A]*(Esoma-p[V_K])*(tmp*x[idx+10]+(1.0-tmp)*x[idx+11]);
  Ih= x[idx+12]*p[g_h]*(Esoma-p[V_h]);
  Il= p[g_leaks]*(Esoma-p[V_leak]);

  os << ICa << " ";
  os << IoCa << " ";
  os << IA << " ";
  os << Ih << " ";
  os << Il << endl;
}

void LPRneuron::derivative(double *x, double *dx)
{
  using namespace LPR5;
  static double a, b, hinf, kh, minf, km;
  static double INa, Id, ICa, IoCa, IA, Ih, Il, IM, IVV;
  static double Isyn;
  
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }

  // differential eqn for the axon membrane potential
  INa= pw3(x[idx+1])*x[idx+2]*p[g_Na]*(Eaxon-p[V_Na]);
  Id= pw4(x[idx+8])*p[g_Kd]*(Eaxon-p[V_K]);
  IM= p[g_M]*x[idx+14]*(Eaxon-p[V_M]);
  Il= p[g_leak]*(Eaxon-p[V_leak]);
  IVV= p[g_VV]*(Esoma-Eaxon);
  
  dx[idx]= (-INa - Id - IM - Il + IVV)/p[C_a]; 

  // soma compartment
  a= exp(Esoma/p[RTF]);
  ICa= (x[idx+3]*x[idx+4]*p[g_CaT]+x[idx+5]*p[g_CaS])
    *p[P_Ca]*Esoma*(x[idx+13]*a-p[Ca_out])/(a-1.0);
  
  IoCa= x[idx+6]*x[idx+7]*p[g_KCa]*(Esoma-p[V_K]);

  a= efunc(Esoma,p[V_aA],p[s_aA]);
  IA= pw3(x[idx+9])*p[g_A]*(Esoma-p[V_K])*(a*x[idx+10]+(1.0-a)*x[idx+11]);

  Ih= x[idx+12]*p[g_h]*(Esoma-p[V_h]);

  Il= p[g_leaks]*(Esoma-p[V_leak]);

  dx[idx+15]= (-ICa - IoCa - IA - Ih - Il - IVV +
	       (p[I_DC]+Isyn)*p[I_scale])/p[C_s];

  // differential eqn for mNa
  a= (3.5+0.1*x[idx]) / (1.0-exp(-3.5-0.1*x[idx]));
  b= 4.0*exp(-(x[idx]+60.0)/18.0);
  dx[idx+1]= (a*(1.0-x[idx+1])-b*x[idx+1])*p[k_slow]*p[k_fast];
  // differential eqn for hNa
  a= 0.07*exp(-x[idx]/20.0-3.0);   
  b= 1.0 / (exp(-3.0-0.1*x[idx])+1.0);
  dx[idx+2]= (a*(1.0-x[idx+2])-b*x[idx+2])*p[k_slow]*p[k_fast];

  // differential eqn for md
  a= (-0.5-0.01*x[idx]) / (exp(-5.0-0.1*x[idx])-1.0); 
  b= 0.125*exp(-(x[idx]+60.0)/80.0);
  dx[idx+8]= (a*(1.0-x[idx+8])-b*x[idx+8])*p[k_slow]*p[k_fast];

  // differential eqn for mCaT
  minf= efunc(Esoma,p[V_mCaT],p[s_mCaT]);
  km= p[k_mCaT]*p[k_slow];
  dx[idx+3]= (minf-x[idx+3])*km;
  // differential eqn for hCaT
  hinf= efunc(Esoma,p[V_hCaT],p[s_hCaT]);
  kh= p[k_hCaT]*efunc(Esoma,p[V_khCaT],p[s_khCaT])*p[k_slow];
  dx[idx+4]= (hinf-x[idx+4])*kh;

  // differential eqn for mCaS
  minf= efunc(Esoma,p[V_mCaS],p[s_mCaS]);
  km= p[k_mCaS]*p[k_slow];
  dx[idx+5]= (minf-x[idx+5])*km;

  // differential eqn for mKCa
  a= p[V_mKCa1]-p[f]*x[idx+13];
  b= p[V_mKca2]-p[f]*x[idx+13];
  minf= efunc(Esoma,a,p[s_mKCa1])
    *efunc(Esoma,b,p[s_mKca2])*(x[idx+13]/(p[c_mKCa]+x[idx+13]));
  km= p[k_mKCa]*p[k_slow];
  dx[idx+6]= (minf-x[idx+6])*km;
  // differential eqn for hKCa
  hinf= p[c_hKCa1]/(p[c_hKCa2]+x[idx+13]);
  kh= p[k_hKCa]*p[k_slow];
  dx[idx+7]= (hinf-x[idx+7])*kh;

  // differential eqn for mA
  minf= efunc(Esoma,p[V_mA],p[s_mA]);
  km= p[k_mA]*p[k_slow];
  dx[idx+9]= (minf-x[idx+9])*km;
  // differential eqn for hA1
  hinf= efunc(Esoma,p[V_hA],p[s_hA]);
  kh= p[k_hA1]*p[k_slow];
  dx[idx+10]= (hinf-x[idx+10])*kh;
  // differential eqn for hA2
  kh= p[k_hA2]*efunc(Esoma,p[V_khA2],p[s_khA2])*p[k_slow];
  dx[idx+11]= (hinf-x[idx+11])*kh;
    
  // differential eqn for mh
  minf= efunc(Esoma,p[V_mh],p[s_mh]); 
  km= p[k_mh]/efunc(Esoma,p[V_kmh],p[s_kmh])*p[k_slow];
  dx[idx+12]= (minf-x[idx+12])*km;
     
  // differential eqn for Ca concentration 
  dx[idx+13]= -p[c_ICa]*ICa - p[k_Ca]*(x[idx+13]-p[Ca_0]);

  // differential equn for M current activation var
  minf= efunc(Eaxon,p[V_mM],p[s_mM]);
  km= p[k_mM]*efunc(Eaxon,p[V_kmM],p[s_kmM])*p[k_slow];
  dx[idx+14]= (minf-x[idx+14])*km;
}

#undef efunc
#undef Eaxon
#undef Esoma

#endif
