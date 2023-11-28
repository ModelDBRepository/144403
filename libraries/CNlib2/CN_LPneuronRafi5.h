//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institute for Nonlinear Dynamics
//            University of California San Diego
//            La Jolla, CA 92093-0402
//
// email to:  tnowotny@ucsd.edu
//
// initial version: 2002-01-24
//
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// adapted from LPneuronJorge with
// a) own fit for Ih
// b) GHK formalism for ICa + own fit with it
//--------------------------------------------------------------------------

#ifndef LPRNEURON_H
#define LPRNEURON_H

#include "CN_neuron.h"
#include <cmath>

#define LPR_PNO 71
#define LPR_IVARNO 16

//--------------------------------------------------------------------------
// parameters of the LPR neuron (typical set, can be modified at any time
// later on. Modifications take effect immediately in the next evaluation)

namespace LPR5
{
  enum {
    g_Na,
    V_Na,
    g_CaT,
    g_CaS,
    g_KCa,
    V_K,
    g_Kd,
    g_A,
    g_h,
    V_h,
    g_leak,
    V_leak,
    C_a,
    c_ICa,
    k_Ca,
    k_slow,
    V_shift,
    I_DC,
    k_fast,
    g_VV,
    C_s,
    g_leaks,
    RTF,
    k_mCaT,
    k_hCaT,
    k_mCaS,
    V_mCaT,
    V_hCaT,
    V_mCaS,
    s_mCaT,
    s_hCaT,
    s_mCaS,
    P_Ca,
    Ca_out,
    V_khCaT,
    s_khCaT,
    k_mKCa,
    k_hKCa,
    V_mKCa1,
    V_mKca2,
    s_mKCa1,
    s_mKca2,
    f,
    c_mKCa,
    c_hKCa1,
    c_hKCa2,
    Ca_0,
    k_mA,
    k_hA1,
    k_hA2,
    V_mA,
    V_hA,
    V_khA2,
    V_aA,
    s_mA,
    s_hA,
    s_khA2,
    s_aA,
    k_mh,
    V_mh,
    V_kmh,
    s_mh,
    s_kmh,
    I_scale,
    V_M,
    g_M,
    k_mM,
    V_mM,
    s_mM,
    V_kmM,
    s_kmM
  };
}

double stdLPR_p[LPR_PNO]= {
  3.7,     // 0 - g_Na
  50,     // 1 - V_Na
  0.331665957,     // 2 - g_CaT
  0.1951207285,     // 3 - g_CaS
  7.64933323,     // 4 - g_KCa
  -72,     // 5 - V_K
  0.3438103179,     // 6 - g_Kd
  1.166565404,     // 7 - g_A
  0.01162925953,     // 8 - g_h
  -20,     // 9 - V_h
  0.001739191812,     // 10 - g_leak
  -50,     // 11 - V_leak
  0.00386062426,     // 12 - C_a
  0.03835925304,     // 13 - c_ICa
  0.001870585507,     // 14 - k_Ca
  1.5382787,     // 15 - k_slow
  -6.125252573,     // 16 - V_shift
  -8.278437745,     // 17 - I_DC
  0.2674443424,     // 18 - k_fast
  0.01051638039,     // 19 - g_VV
  0.2512112285,     // 20 - C_s
  0.001990003389,     // 21 - g_leaks
  11.4866,     // 22 - RTF
  0.045,     // 23 - k_mCaT
  0.02,     // 24 - k_hCaT
  0.0608644,     // 25 - k_mCaS
  15,     // 26 - V_mCaT
  -40,     // 27 - V_hCaT
  29.1336,     // 28 - V_mCaS
  -9.8,     // 29 - s_mCaT
  3.2,     // 30 - s_hCaT
  -4.43083,     // 31 - s_mCaS
  0.00155982,     // 32 - P_Ca
  15123.5,     // 33 - Ca_out
  -15,     // 34 - V_khCaT
  -10,     // 35 - s_khCaT
  0.6,     // 36 - k_mKCa
  0.035,     // 37 - k_hKCa
  0,     // 38 - V_mKCa1
  -16,     // 39 - V_mKca2
  -23,     // 40 - s_mKCa1
  -5,     // 41 - s_mKca2
  0.6,     // 42 - f
  2.5,     // 43 - c_mKCa
  0.7,     // 44 - c_hKCa1
  0.6,     // 45 - c_hKCa2
  0.001728873164,     // 46 - Ca_0
  0.14,     // 47 - k_mA
  0.05,     // 48 - k_hA1
  0.002519870695,     // 49 - k_hA2
  -12,     // 50 - V_mA
  -62,     // 51 - V_hA
  -40,     // 52 - V_khA2
  7,     // 53 - V_aA
  -26,     // 54 - s_mA
  6,     // 55 - s_hA
  -12,     // 56 - s_khA2
  -15,     // 57 - s_aA
  5.5e-05,     // 58 - k_mh
  -70,     // 59 - V_mh
  -110,     // 60 - V_kmh
  8,     // 61 - s_mh
  -21.6,     // 62 - s_kmh
  0.01967580907,     // 63 - I_scale
  -80,     // 64 - V_M
  0.493711703,     // 65 - g_M
  0.002424418911,     // 66 - k_mM
  -38.33950896,     // 67 - V_mM
  -3.973413716,     // 68 - s_mM
  -50.11214269,     // 69 - V_kmM
  15.17576767,     // 70 - s_kmM
};

double *LPR_p= stdLPR_p;

const char *LPR_p_text[LPR_PNO]= {
  "g_Na",
  "V_Na",
  "g_CaT",
  "g_CaS",
  "g_KCa",
  "V_K",
  "g_Kd",
  "g_A",
  "g_h",
  "V_h",
  "g_leak",
  "V_leak",
  "C_a",
  "c_ICa",
  "k_Ca",
  "k_slow",
  "V_shift",
  "I_DC",
  "k_fast",
  "g_VV",
  "C_s",
  "g_leaks",
  "RTF",
  "k_mCaT",
  "k_hCaT",
  "k_mCaS",
  "V_mCaT",
  "V_hCaT",
  "V_mCaS",
  "s_mCaT",
  "s_hCaT",
  "s_mCaS",
  "P_Ca",
  "Ca_out",
  "V_khCaT",
  "s_khCaT",
  "k_mKCa",
  "k_hKCa",
  "V_mKCa1",
  "V_mKca2",
  "s_mKCa1",
  "s_mKca2",
  "f",
  "c_mKCa",
  "c_hKCa1",
  "c_hKCa2",
  "Ca_0",
  "k_mA",
  "k_hA1",
  "k_hA2",
  "V_mA",
  "V_hA",
  "V_khA2",
  "V_aA",
  "s_mA",
  "s_hA",
  "s_khA2",
  "s_aA",
  "k_mh",
  "V_mh",
  "V_kmh",
  "s_mh",
  "s_kmh",
  "I_scale",
  "V_M",
  "g_M",
  "k_mM",
  "V_mM",
  "s_mM",
  "V_kmM",
  "s_kmM"
};


//----------------------------------------------------------------------
// the following initial values reflect the steady state with no input

double stdLPR_INIVARS[LPR_IVARNO]= {
  -56.0,                       // 0 - membrane potential 
  0.0,                         // 1 - Na channel activation mNa
  0.2,                         // 2 - Na channel unblocking hNa
  0.0,                         // 3 - CaT channel activation mCaT
  1.0,                         // 4 - CaT channel unblocking hCaT
  0.0,                         // 5 - CaS channel activation mCaS
  0.0,                         // 6 - oCa channel activation moCa
  1.0,                         // 7 - oCa channel unblocking hoCa
  0.0,                         // 8 - d channel activation md
  0.0,                         // 9 - A channel activation mA
  1.0,                         // 10 - A channel unblocking hA1
  1.0,                         // 11 - A channel unblocking hA2
  0.0,                         // 12 - Ih channel activation
  0.04,                        // 13 - Ca concentration
  0.01,                        // 14 - M current activation var
  -56.0,                       // 15 - membrane potential soma
};

double *LPR_INIVARS= stdLPR_INIVARS;

const char *LPR_INIVARSTEXT[LPR_IVARNO]= {
  "0 - membrane potential",
  "1 - Na channel activation mNa",
  "2 - Na channel unblocking hNa",
  "3 - CaT channel activation mCaT",
  "4 - CaT channel unblocking hCaT",
  "5 - CaS channel activation mCaS",
  "6 - oCa channel activation moCa",
  "7 - oCa channel unblocking hoCa",
  "8 - d channel activation md",
  "9 - A channel activation mA",
  "10 - A channel unblocking hA1",
  "11 - A channel unblocking hA2",
  "12 - Ih channel activation",
  "13 - Ca concentration",
  "14 - M current activation var",
  "15 - membrane potential soma"
};


// for fitting purposes: type == 0 means additive changes
//                       type == 1 means multiplicative changes

double LPR_p_type[LPR_PNO]= {
  1,      // 0 - g_Na
  0,      // 1 - V_Na
  1,      // 2 - g_CaT
  1,      // 3 - g_CaS
  1,      // 4 - g_KCa
  0,      // 5 - V_K
  1,      // 6 - g_Kd
  1,      // 7 - g_A
  1,      // 8 - g_h
  0,      // 9 - V_h
  1,      // 10 - g_leak
  0,      // 11 - V_leak
  1,      // 12 - C_a
  1,      // 13 - c_ICa
  1,      // 14 - k_Ca
  1,      // 15 - k_slow
  0,      // 16 - V_shift
  0,      // 17 - I_DC
  1,      // 18 - k_fast
  1,      // 19 - g_VV
  1,      // 20 - C_s
  1,      // 21 - g_leaks
  1,      // 22 - RTF
  1,      // 23 - k_mCaT
  1,      // 24 - k_hCaT
  1,      // 25 - k_mCaS
  0,      // 26 - V_mCaT
  0,      // 27 - V_hCaT
  0,      // 28 - V_mCaS
  1,      // 29 - s_mCaT
  1,      // 30 - s_hCaT
  1,      // 31 - s_mCaS
  1,      // 32 - P_Ca
  1,      // 33 - Ca_out
  0,      // 34 - V_khCaT
  1,      // 35 - s_khCaT
  1,      // 36 - k_mKCa
  1,      // 37 - k_hKCa
  0,      // 38 - V_mKCa1
  0,      // 39 - V_mKca2
  1,      // 40 - s_mKCa1
  1,      // 41 - s_mKca2
  1,      // 42 - f
  1,      // 43 - c_mKCa
  1,      // 44 - c_hKCa1
  1,      // 45 - c_hKCa2
  1,      // 46 - Ca_0
  1,      // 47 - k_mA
  1,      // 48 - k_hA1
  1,      // 49 - k_hA2
  0,      // 50 - V_mA
  0,      // 51 - V_hA
  0,      // 52 - V_khA2
  0,      // 53 - V_aA
  1,      // 54 - s_mA
  1,      // 55 - s_hA
  1,      // 56 - s_khA2
  1,      // 57 - s_aA
  1,      // 58 - k_mh
  0,      // 59 - V_mh
  0,      // 60 - V_kmh
  1,      // 61 - s_mh
  1,      // 62 - s_kmh
  1,      // 63 - I_scale
  0,      // 64 - V_M
  1,      // 65 - g_M
  1,      // 66 - k_mM
  0,      // 67 - V_mM
  1,      // 68 - s_mM
  0,      // 69 - V_kmM
  1       // 70 - s_kmM
};

// for fitting purposes ... the LPR_p_sens(itivity) is basically the inverse
// of the slope of the deviation of the neuron model perturbed in a given
// parameter compared to the model with the original parameter set.
// see "test_params" and "calc_sensitivity" for more details

double LPR_p_sens[LPR_PNO]= {
  5e-10,     // 0 - g_Na
  1e-09,     // 1 - V_Na
  1.56905298e-10,     // 2 - g_CaT
  1.725958278e-10,     // 3 - g_CaS
  2.95245e-10,     // 4 - g_KCa
  1e-09,     // 5 - V_K
  5e-10,     // 6 - g_Kd
  2.15233605e-10,     // 7 - g_A
  2.796052411e-10,     // 8 - g_h
  1e-09,     // 9 - V_h
  5.31441e-10,     // 10 - g_leak
  1e-09,     // 11 - V_leak
  5.845851e-10,     // 12 - C_a
  5.9049e-10,     // 13 - c_ICa
  5.9049e-10,     // 14 - k_Ca
  1.398026206e-10,     // 15 - k_slow
  4.73513931e-10,     // 16 - V_shift
  1.412147682e-09,     // 17 - I_DC
  2.3914845e-10,     // 18 - k_fast
  2.541865828e-10,     // 19 - g_VV
  3.486784401e-10,     // 20 - C_s
  2.51644717e-10,     // 21 - g_leaks
  1e-09,     // 22 - RTF
  1e-09,     // 23 - k_mCaT
  1e-09,     // 24 - k_hCaT
  1e-09,     // 25 - k_mCaS
  1e-09,     // 26 - V_mCaT
  1e-09,     // 27 - V_hCaT
  1e-09,     // 28 - V_mCaS
  1e-09,     // 29 - s_mCaT
  1e-09,     // 30 - s_hCaT
  1e-10,     // 31 - s_mCaS
  1e-09,     // 32 - P_Ca
  1e-09,     // 33 - Ca_out
  5e-10,     // 34 - V_khCaT
  5e-10,     // 35 - s_khCaT
  1e-09,     // 36 - k_mKCa
  5e-10,     // 37 - k_hKCa
  1e-09,     // 38 - V_mKCa1
  1e-09,     // 39 - V_mKca2
  5e-09,     // 40 - s_mKCa1
  1e-09,     // 41 - s_mKca2
  1e-09,     // 42 - f
  5e-10,     // 43 - c_mKCa
  1e-09,     // 44 - c_hKCa1
  1e-09,     // 45 - c_hKCa2
  1e-09,     // 46 - Ca_0
  1e-09,     // 47 - k_mA
  2.63063295e-09,     // 48 - k_hA1
  1e-09,     // 49 - k_hA2
  1e-09,     // 50 - V_mA
  1e-09,     // 51 - V_hA
  1e-09,     // 52 - V_khA2
  1e-09,     // 53 - V_aA
  1e-09,     // 54 - s_mA
  1e-09,     // 55 - s_hA
  5e-10,     // 56 - s_khA2
  1e-09,     // 57 - s_aA
  1e-09,     // 58 - k_mh
  1e-09,     // 59 - V_mh
  5e-09,     // 60 - V_kmh
  1e-09,     // 61 - s_mh
  1e-09,     // 62 - s_kmh
  1e-09,     // 63 - I_scale
  1e-09,     // 64 - V_M
  1.66771817e-10,     // 65 - g_M
  1e-09,     // 66 - k_mM
  2.824295365e-10,     // 67 - V_mM
  3.486784401e-10,     // 68 - s_mM
  3.486784401e-10,     // 69 - V_kmM
  1.500946353e-10,     // 70 - s_kmM
};


// the LPR neuron class itself

class LPRneuron: public neuron
{
 public:
  LPRneuron(int, double *);
  LPRneuron(int, vector<int>, double *);
  virtual ~LPRneuron() { }
  inline virtual double E(double *);
  virtual void currents(ostream &, double *);
  virtual void derivative(double *, double *);
};


#endif
