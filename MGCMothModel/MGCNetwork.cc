/*--------------------------------------------------------------------------
 Author: Christopher L Buckley and Thomas Nowotny

 Institute: Centre for Computational Neuroscience and Robotics
 University of Sussex
 FcondNetworkmer, Brighton BN1 9QJ, UK

 email to:  c.l.buckley@sussex.ac.uk

 version: 2009-11-02
 --------------------------------------------------------------------------*/

using namespace std;

#include<cstdlib>
#include <sstream>
#include "CN_DCInput.h"
#include "CN_Rallsynapse.h"
#include "CN_absynapse.h"
#include "CN_rk65n.h"
#include "tnt.h"
#include <complex>
#include "jama_eig.h"
#include <vector>
#include <sys/time.h>
#include "gauss.h"

stdRG RG;
randomGauss RGaus;
#include "CN_Emptysynapse.h"
#include "CN_NeuronModel.cc"
#include "CN_absynapse.cc"
#include "CN_DCInput.cc"
#include "CN_Emptysynapse.cc"
#include "CN_rk65n.cc"
#include <iostream>
#include <fstream>

//Flags to set
bool doCond = false;
bool doPatterned =  false;

int stimNum = 10000;
bool doRate = true;
bool backwardScale = true;
bool Manual = true;
bool CreateFile = true;
bool CaAdapt = false;

bool doneFileCreate = true;

//the parameters
double gm = 1.0;
double tauz = 0.02;
double betaGlobal = 0.005;
double Ic = 0.04132;

double minFreq = 0.008;
double maxFreq = 0.020;

double NoiseMag = 0.000;

//this are the parameter to gradually scale the weights over scaleTime ms in bins of length scaleBin
double scaleBin = 1;
double scaleTime = 5000;
double scaleStart2 = 6000;
double scaleTime2 = 8000;


double percentCritical = 0.85;

double INPUTMAG = 0.1;
int THETIME = 20000;
//int THETIME = 8000;
int IMPULSESTART = 8000;
//int IMPULSESTART = 2000;
double TURNOFFADAPT;
int IMPULSEDUR = 4000;
double OutTime;
//Manualy set

int NoFast = 0;
int NoInput =0;

int NoPNs = 20;
int NoORNs = 1;
int NoLNs = 20;



double ORN2LNWeight = 0.0006;
	double ORN2PNWeight = 0.00055;
	double LN2PNWeight = -30;

double LN2PNProb = 0.5;
double Pcon = 0.75;




double gSS = 1;
double gUU = 1;
double gSU = 1;// the stimulated   nodes forward connections
double gUS = 1;

int NoTotal = NoPNs + NoLNs;


char *globalName;

#define NEUROSYNADAPT_PNO 19
#define NEUROSYNADAPT_IVARNO 8
#define NEUROSYN_PNO 19

Array1D<double> SeqArray;
double Seq;
double SeqMax;
bool SlowerAdapt = false;

//Flag to output data: Always true for single network run
bool OutDat = true;
bool Gradual;
//Linear fit of how the rest potential depends on the input current
double a;
double b;
double c;

//Linear fit to the F-I curve
double M;
double C;

//Ermentraoyt fit
double A;
double B;
double A2;
Array1D<double> vrestArray;

double PARAMS_LN[NEUROSYNADAPT_PNO] = { 7.15, // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
		50.0, // 1 - ENa: Na equi potential in mV
		1.43, // 2 - gK: K conductance in 1/(mOhms * cm^2)
		-95.0, // 3 - EK: K equi potential in mV
		0.021, // 4 - gl: leak conductance in 1/(mOhms * cm^2)
		-55.0, // 5 - El: leak equi potential in mV
		0.00572, // 6 - gKl: potassium leakage conductivity
		-95.0, // 7 - EKl: potassium leakage equi pot in mV
		65.0, // 8 - V0: ~ total equi potential (?)
		0.143, // 9 - Cmem: membr. capacity density in muF/cm^2
		gm,//0.715, // 10 - gM: conductance of the M current
		0.0, // 11- IDC: baseline offset current
		-95, //  12 inEsyn %reversal potential (-95 = inhibitory)
		-20, // 13  inEpre %threshold for pre synaptic spike detection
		2, //14 inasyn
		betaGlobal, // 15 inbsyn
		1, // 16 inrtime
		0, //17 noise
		0 //Equilibrium Synspe value
		};

double PARAMS_LNFAST[NEUROSYNADAPT_PNO] = { 7.15, // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
		50.0, // 1 - ENa: Na equi potential in mV
		1.43, // 2 - gK: K conductance in 1/(mOhms * cm^2)
		-95.0, // 3 - EK: K equi potential in mV
		0.021, // 4 - gl: leak conductance in 1/(mOhms * cm^2)
		-55.0, // 5 - El: leak equi potential in mV
		0.00572, // 6 - gKl: potassium leakage conductivity
		-95.0, // 7 - EKl: potassium leakage equi pot in mV
		65.0, // 8 - V0: ~ total equi potential (?)
		0.143, // 9 - Cmem: membr. capacity density in muF/cm^2
		gm,//0.715, // 10 - gM: conductance of the M current
		0.0, // 11- IDC: baseline offset current
		-95, //  12 inEsyn %reversal potential (-95 = inhibitory)
		-20, // 13  inEpre %threshold for pre synaptic spike detection
		2, //14 inasyn
		0.1, // 15 inbsyn
		1, // 16 inrtime
		0, //17 noise
		0 //Equilibrium Synspe value
		};

double PARAMS_PN[NEUROSYNADAPT_PNO] = { 7.15, // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
		50.0, // 1 - ENa: Na equi potential in mV
		1.43, // 2 - gK: K conductance in 1/(mOhms * cm^2)
		-95.0, // 3 - EK: K equi potential in mV
		0.021, // 4 - gl: leak conductance in 1/(mOhms * cm^2)
		-55.0, // 5 - El: leak equi potential in mV
		0.00572, // 6 - gKl: potassium leakage conductivity
		-95.0, // 7 - EKl: potassium leakage equi pot in mV
		65.0, // 8 - V0: ~ total equi potential (?)
		0.143, // 9 - Cmem: membr. capacity density in muF/cm^2
		gm,//0.715, // 10 - gM: conductance of the M current
		0.0, // 11- IDC: baseline offset current
		0, //  12 inEsyn %reversal potential (-95 = inhibitory)
		-20, // 13  inEpre %threshold for pre synaptic spike detection
		2, //14 inasyn
		betaGlobal, // 15 inbsyn
		1, // 16 inrtime
		0, //17 noise
		0 //Equilibrium Synspe value
		};

double PARAMS_ORN[NEUROSYNADAPT_PNO] = { 2, // 0 - spike time of mf_poisson inputneuron
		10.0, // 1 - refractory period + spike time
		-60.0, // 2 - input neuron resting potential
		50.0, // 3 - input neuron potential when firing
		0.0, // 4 - firing rate Lambda [1/ms]=[10^3 Hz]
		2, //5 inasyn
		0.01, // 6 inbsyn
		1, // 7 inrtime
		-55.0, // 8
		0, // 9
		0,//0.715, // 10
		0.0, // 11-
		0.0, //  12
		0, // 13
		0, //14
		0, // 15
		0, // 16
		0, //17
		0 //
		};

const char *NEUROSYN_p_text[NEUROSYNADAPT_PNO] = {
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
		"12 - gKl: potassium leakage conductivity",
		"13 - EKl: potassium leakage equi pot in mV",
		"14 - V0: ~ total equi potential (?)",
		"15 - Cmem: membr. capacity density in muF/cm^2",
		"16 - gM: conductance of the M current",
		"17 - IDC: baseline offset current"
			"18 - eq valuet" };

double NEUROSYNADAPT_INIVARS[NEUROSYNADAPT_IVARNO] = { -60.0, // 0 - membrane potential E
		0.0529324, // 1 - prob. for Na channel activation m
		0.3176767, // 2 - prob. for not Na channel blocking h
		0.5961207, // 3 - prob. for K channel activation n
		0.1, // 4 - M current activation
		0, // 5 - synapse activation
		0.0, // 6 - calcium current
		0.0 // 7 - that current
		};

const char *NEUROSYN_INIVARSTEXT[NEUROSYNADAPT_IVARNO] = {
		"0 - membrane potential E", "1 - prob. for Na channel activation m",
		"2 - prob. for not Na channel blocking h",
		"3 - prob. for K channel activation n", "4 - M current activation",
		"5 - synapse activation", "6 - calcium" };

#include "CN_NeuroSynAdapt.h"
#include "CN_NeuroSynRate.h"
#include "CN_NeuroSynAdapt.cc"
#include "CN_NeuroSynRate.cc"
#include "CN_PoissonRateNeuron.h"
#include "CN_PoissonRateNeuron.cc"
#include "CN_Poissonneuron.h"
#include "CN_Poissonneuron.cc"

double rk65_MINDT = 0.05;
double rk65_eps = 1e-12;
double rk65_relEps = 1e-9;
double rk65_absEps = 1e-16;
double mindt = 1e-6;

#define EXTRAWRITE cerr << x[0] <<" "; n.currents(cerr,x);

#include "fullMGC.h"
#include "fullMGCRate.h"

int main(int argc, char *argv[]) {


	//data output file
	stringstream rname, cname;
	ofstream rateEnd, memEnd;
	rateEnd.precision(10);
	memEnd.precision(10);

	char thername[80];
	char thecname[80];

	if (tauz == 0.05 && gm == 1.0) {
		//0.05{
		A = 0.0848;
		B = 13.1848;
		a = -0.5678;
		b = 5.8296;
		c = -67.4830;

		M = 0.0347;
		C = 0.0114;
	} else if (tauz == 0.02 && gm == 1.0) {
		//0.02
		A = 0.1224;
		B = 21.1075;

		a = -0.5967;
		b = 5.2418;
		c = -65.0945;

		M = 0.0353;
		C = 0.0043;
	} else if (tauz == 0.01 && gm == 1.0) {
		//0.01
		A = 0.1672;
		B = 24.0765;
		a = -0.4664;
		b = 4.4276;
		c = -63.5274;

		M = 0.0359;
		C = 0.0015;
	} else if (tauz == 0.01 && gm == 2.0) {
		A = 0.0632;
		B = 39.8488;
		a = -0.4664;
		b = 4.4276;
		c = -63.5274;

		M = 0.0185;
		C = 0.0023;
	}

	A2 = pow(A, 2);
	int seeder = 0;

	if (Manual) {
		if (argc != 23) {
			cerr
					<< "usage: filename percentCritical NoLNs NoPNs Pcon NoORNs  gSS gUS gSU gUU minFreq maxFreq  backwardScale docond ORN2LNWeight doneFileCreate LN2PNWeight ORN2PNWeight IMPULSEDUR stimNum IMPULSESTART seed"
					<< endl;
			exit(1);
		}

		cerr << "call was: ";
		for (int i = 0; i < argc; i++) {
			cerr << argv[i] << " ";
		}
		cerr << endl;

		globalName = argv[1];

		percentCritical = atof(argv[2]);
		cout << "percentCritical: " << percentCritical << " ";

		NoLNs = atoi(argv[3]);
		cout << "NoLNs: " << NoLNs << " ";

		NoPNs = atoi(argv[4]);
				cout << "NoPNs: " << NoPNs << " ";

		Pcon = atof(argv[5]);
		stringstream name;
		cout << "Pcon: " << Pcon << " ";

		NoORNs = atoi(argv[6]);
		cout << "NoORNs: " << NoORNs << " ";

		gSS = (atof(argv[7]));
		cout << "gSS: " << gSS << " ";

		gUS = (atof(argv[8]));
		cout << "gUS: " << gUS << " ";

		gSU = (atof(argv[9]));
		cout << "gSU: " << gSU << " ";

		gUU = (atof(argv[10]));
		cout << "gUU: " << gUU << " ";

		minFreq = (atof(argv[11]));
		cout << "minFreq: " << minFreq << " ";

		maxFreq = (atof(argv[12]));
		cout << "maxFreq: " << maxFreq << " ";

		backwardScale = bool(atoi(argv[13]));
		cout << "backwardScale: " << backwardScale << " ";

		doCond = atoi(argv[14]);
		cout << "doCond: " << doCond << " ";

		ORN2LNWeight = (atof(argv[15]));
		cout << "ORN2LNWeight: " << ORN2LNWeight << " ";

		doneFileCreate = atoi(argv[16]);
		cout << "doneFileCreate: " << doneFileCreate << " ";

		LN2PNWeight = atof(argv[17]);
					cout << "LN2PNWeight: " << LN2PNWeight << " ";

		ORN2PNWeight = atof(argv[18]);
				cout << "ORN2PNWeight: " << ORN2PNWeight << " ";


		IMPULSEDUR = atoi(argv[19]);
				cout << "IMPULSEDUR: " << IMPULSEDUR << " ";



				stimNum = atoi(argv[20]);
						cout << "stimNum: " << stimNum << " ";



						IMPULSESTART = atoi(argv[21]);
								cout << "IMPULSESTART: " << IMPULSESTART << " ";


		seeder = atoi(argv[22]);
		cout << "seed: " << seeder << " ";


		RG.seed(seeder);

		rname << argv[1] << "rate.dat";
		cname << argv[1] << "cond.dat";

	     NoTotal = NoPNs + NoLNs;

	} else {
		seeder = int(RG.n() * 100);
		//seeder=17;
		RG.seed(seeder);
		rname /*<< int(log(INPUTMAG)*100)*/<< "rate.dat";
		cname /* << int(log(INPUTMAG)*100)*/<< "cond.dat";
	}

	LN2PNWeight = LN2PNWeight/(LN2PNProb*NoLNs);
	ORN2PNWeight = ORN2PNWeight/NoORNs;
	if(doCond)
		 Gradual = true;
	else
		 Gradual = false;

	if(Gradual)
	{

		THETIME = THETIME+scaleStart2+scaleTime2+2000;
		IMPULSESTART = IMPULSESTART+scaleStart2+scaleTime2+2000;
		 TURNOFFADAPT = IMPULSESTART-5000;
	}

	OutTime = IMPULSESTART-4000;
	rname >> thername;
	cname >> thecname;
	rateEnd.open(thername);
	memEnd.open(thecname);

	//Seed


	//Set equilibrium postions
	Array1D<double> SeqArrayInit(NoTotal, 0.0);
	SeqArray = SeqArrayInit;

	cout << " The synaptic eqm values are" << endl;
	for (int i = 0; i < NoTotal; i++) {

		if (i < NoFast) {
			double alpha = PARAMS_LNFAST[14];
			double beta = PARAMS_LNFAST[15];
			double tr = PARAMS_LNFAST[16];
			Seq = alpha * tr * minFreq / beta;
			SeqMax = alpha * tr * maxFreq / beta;
		}

		if (i >= NoFast && i < NoLNs) {
			double alpha = PARAMS_LN[14];
			double beta = PARAMS_LN[15];
			double tr = PARAMS_LN[16];
			Seq = alpha * tr * minFreq / beta;
			SeqMax = alpha * tr * maxFreq / beta;
		}

		if (i>=NoLNs)
			{
				double alpha = PARAMS_LN[14];
				double beta = PARAMS_LN[15];
				double tr = PARAMS_LN[16];
				Seq = alpha * tr * minFreq / beta;
						SeqMax = alpha * tr * maxFreq / beta;

			}

			SeqArray[i] = Seq + (SeqMax - Seq) * RG.n();
			cout << SeqArray[i] << " ";
		}
		cout << endl;

		//the networks

		fullMGCRate NSRateNet;
		fullMGC NSNet;

		Array2D<double> weightsRate(NoTotal, NoTotal, 0.0);
		Array2D<double> weightsCond(NoTotal, NoTotal, 0.0);
		//set weights with eual incomign number


		int NoCon = int(NoLNs * Pcon);
		for (int i = 0; i < NoLNs; i++) {

			weightsRate[i][i] = 0.0;
			int counter = 0;
			Array1D<double> holdInt(NoCon, 99999.0);

			while (counter < NoCon) {
				int tryInt = int(NoLNs * RG.n());

				bool flag = true;
				for (int j = 0; j < NoCon; j++) {
					if (holdInt[j] == tryInt)
						flag = false;
				}

				if (tryInt != i && flag) {
					weightsRate[tryInt][i] = -1.0;
					holdInt[counter] = tryInt;
					counter++;
				}
			}
		}

		for (int i = 0; i < NoLNs; i++) {
			for (int j = 0; j < NoLNs; j++) {

				//stimulated node interconnections
				if (j < NoORNs && i < NoORNs)
					weightsRate[j][i] = weightsRate[j][i] * gSS;

				// the stimulated   nodes forward connections
				if (j < NoORNs && i >= NoORNs)
					weightsRate[j][i] = weightsRate[j][i] * gSU;

				//unstimulated node interconnections
				if (j >= NoORNs && i >= NoORNs)
					weightsRate[j][i] = weightsRate[j][i] * gUU;

				// the ustimulated   nodes backward connections
				if (j >= NoORNs && i <= NoORNs)
					weightsRate[j][i] = weightsRate[j][i] * gUS;

				if (i == j)
					weightsRate[j][i] = 0.0;
			}

		}

		for (int i = NoLNs; i < NoTotal; i++) {
			for (int j = 0; j < NoLNs; j++) {

				if (RG.n() < LN2PNProb)
					weightsRate[j][i] = LN2PNWeight;
			}
		}

		NSRateNet.SetWeights(weightsRate, SeqArray);
		weightsCond = NSRateNet.SetCritical(percentCritical);

		//for turning off the backward connections
		if (!backwardScale) {
			for (int i = 0; i < NoTotal; i++) {
				for (int j = 0; j < NoTotal; j++) {

					if (j < NoORNs && i >= NoORNs)
						weightsCond[j][i] = weightsCond[j][i];
					else
						weightsCond[j][i] = 0.0;

				}
			}
		}

		Array1D<double> biasCond(NoTotal, 0.0);
		NSRateNet.SetWeights(weightsCond, SeqArray);
		biasCond = NSRateNet.SetCC();

		NSNet.SetWeights(weightsCond, SeqArray);
		NSNet.SetCC(biasCond);

		NSRateNet.PrintWeights();
		NSRateNet.PrintThetaValues();
		NSRateNet.PrintMaxEig();

		//trun on input at same time
		Array1D<int> inputStart(NoORNs, 0);
		for (int i = 0; i < NoORNs; i++)
			inputStart[i] = 0;

		Array1D<double> endPointsR(2 * NoTotal, 0.0);

		Array1D<double> endPointsM(2 * NoTotal, 0.0);

		if (doRate) {
			NSRateNet.init();
			endPointsR = NSRateNet.run(THETIME, INPUTMAG, inputStart);
		}
		rateEnd << INPUTMAG << " ";

		for (int i = 0; i < 2 * NoTotal; i++)
			rateEnd << endPointsR[i] << " ";

		rateEnd << "\n";
		rateEnd.close();

		if (doCond) {
			NSNet.init();
			endPointsM = NSNet.run(THETIME, INPUTMAG, inputStart);
		}

		memEnd << INPUTMAG << " ";

		for (int i = 0; i < 2 * NoTotal; i++)
			memEnd << endPointsM[i] << " ";

		memEnd << "\n";
		memEnd.close();

		cout << endl << "The random seed:" << seeder;

		return 0;

	}
