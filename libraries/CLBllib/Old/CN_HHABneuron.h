//--------------------------------------------------------------------------
// Author: Christopher L. Buckley
//
// Institute: CCNR
//			  Informatics
//            University of Sussex
//
// email to:  c.l.buckley@sussex.ac.uk
//
// initial version: 2009-12-15
//
//--------------------------------------------------------------------------


#ifndef CN_HHABNEURON_H
#define CN_HHABNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define HHAB_IVARNO 1
#define HHAB_NPNO 12
#define HHAB_SPNO 7
#define HHAB_PNO 19



double stdHHAB_p[HHAB_PNO]= {
		  0,          // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
		  0,          // 1 - ENa: Na equi potential in mV
		  0,          // 2 - gK: K conductance in 1/(mOhms * cm^2)
		  0,         // 3 - EK: K equi potential in mV
		  0,         // 4 - gl: leak conductance in 1/(mOhms * cm^2)
		  0,         // 5 - El: leak equi potential in mV
		  0,       // 6 - gKl: potassium leakage conductivity
		  0,         // 7 - EKl: potassium leakage equi pot in mV
		  0,          // 8 - V0: ~ total equi potential (?)
		  0,         // 9 - Cmem: membr. capacity density in muF/cm^2
		  0,           // 10 - gM: conductance of the M current
		  0,           // 11- IDC: baseline offset current
		  0, //  ingsyn
		  		 0, //  inEsyn %reversal potential
		  		0, //  inEpre %threshold for pre synaptic spike detection
		  		 0, // inasyn
		  		 0, // inbsyn
		  		0, // inrtime
		  		0 //noise

		};


double *HHAB_p= stdHHAB_p;

const char *HHAB_p_text[HHAB_PNO] = {
		"0 - gNa: Na conductance in 1/(mOhms * cm^2)",
		"1 - ENa: Na equi potential in mV",
		"2 - gK: K conductance in 1/(mOhms * cm^2)",
		"3 - EK: K equi potential in mV",
		"4 - gl: leak conductance in 1/(mOhms * cm^2)",
		"5 - El: leak equi potential in mV",
		"6 - gKl: potassium leakage conductivity",
		"7 - EKl: potassium leakage equi pot in mV",
		"8 - V0: ~ total equi potential (?)",
		"9 - Cmem: membr. capacity density in muF/cm^2",
		"10 - gM: conductance of the M current",
		"11 - IDC: baseline offset current",
		"12 - ingsyn",
		"13 - inEsyn %reversal potential",
		"14 - inEpre %threshold for pre synaptic spike detection",
		"15 - inasyn",
		"16 - inbsyn",
		"17 - inrtime",
		"18 -noise " };


double HHAB_INIVARS[HHAB_IVARNO]= {
  0.1                       // 0 - firing rate
};

const char *HHAB_INIVARSTEXT[HHAB_IVARNO]= {
  "0 - firing rate"
};


// HHAB neuron class itself

class HHABneuron: public neuron
{
 private:
  double Isyn;
 public:
  HHABneuron(int, double *);
  HHABneuron(int, vector<int>, double *);
  ~HHABneuron() { }
  inline virtual double E(double *);
  virtual void set_p(double *, double *);
  virtual void derivative(double *, double *);
  virtual void noise(double *, double *);
};

#endif



