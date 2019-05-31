/*--------------------------------------------------------------------------
 Author: Christopher L Buckley

 Institute: Centre for Computational Neuroscience and Robotics
 University of Sussex
 FcondNetworkmer, Brighton BN1 9QJ, UK

 email to:  c.l.buckley@sussex.ac.uk

 version: 2009-11-02
 --------------------------------------------------------------------------*/

#ifndef NEUROSYNNETWORKRATE_H
#define NEUROSYNNETWORKRATE_H


//w[i][j] source i  and target j

class neuroSynNetworkRate {
public:
	list<neuron *> neurs;
	list<synapse *> syns;
	list<neuron *>::iterator niter;
	list<synapse *>::iterator siter;

	NeuroSynRate **theLNs;
	Emptysynapse **theSynapses;

	DCInput **directLNInput;
	NeuronModel *model;
	rk65n *machine;
	double *x, *xn;
	double dt, dtx;
	int N;

	Array2D<double> mWeights;
	Array1D<double> mSeqArray;
	Array1D<double> vrestMod;
	Array1D<double> mCCBias;
	void GenNetwork(void);

	neuroSynNetworkRate();
	~neuroSynNetworkRate();

	void generateNetwork();
	void generateConnect();
	void init();
	void enable();
	void PrintThetaValues();
	void PrintWeights();
	void PrintMaxEig();

	Array1D<double> run(double tme, double inputCurr,Array1D<int> inputStart);
	double TransferInv(double in, double *params);
	double Transfer(double in, double *params);
	double TransferDeriv(double in, double *params);
	Array2D<double> getJacobian(void);
	void SetWeights(Array2D<double> weights, Array1D<double> SeqArray);
	void ScaleWeights(double factor);
	Array2D<double> SetCritical(double epsilon);
	Array1D<double> SetCC();
	double GetMaxEig();
};

neuroSynNetworkRate::neuroSynNetworkRate() {
	dt = 0.0001;
	dtx = 0.0;
	N = 0.0;

	directLNInput = new DCInput*[NoLNs];
	theLNs = new NeuroSynRate*[NoLNs];


	//Set up network memory and space


	for (int i = 0; i < NoLNs; i++) {
		theLNs[i] = new NeuroSynRate(i, PARAMS_LN);
		neurs.push_back(theLNs[i]);
		directLNInput[i] = new DCInput(theLNs[i], 0);
	}



	//set up connections
	theSynapses = new Emptysynapse*[NoLNs * NoLNs];
	for (int i = 0; i < NoLNs; i++) {
		for (int j = 0; j < NoLNs; j++) {

			if (i != j) {
				theSynapses[j * NoLNs + i] = new Emptysynapse(theLNs[i],
						theLNs[j], 0.0);
			}
		}
	}



	//enable integaror
	enable();
}

neuroSynNetworkRate::~neuroSynNetworkRate() {
	list<neuron *>::iterator i;
	list<synapse *>::iterator j;
	for (i = neurs.begin(); i != neurs.end(); i++) {
		delete *i;
	}
	for (j = syns.begin(); j != syns.end(); j++) {
		delete *j;
	}

	delete[] theLNs;
	delete[] directLNInput;

}

void neuroSynNetworkRate::enable() {
	model = new NeuronModel(&neurs, &syns, N, cerr);
	x = new double[N];
	xn = new double[N];
	machine = new rk65n(N, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);
}

void neuroSynNetworkRate::SetWeights(Array2D<double> weights,
		Array1D<double> SeqArray) {

	mWeights = weights;
	mSeqArray = SeqArray;

	double beta = PARAMS_LN[15];


	//this code alters the wigtsh absed on the rest potetial.

	Array1D<double> setter(NoLNs, 0.0);

	vrestMod = setter;

	double P2[NEUROSYN_PNO];

	for (int i = 0; i < NEUROSYN_PNO; i++)
		P2[i] = PARAMS_LN[i];

	for (int i = 0; i < NoLNs; i++) {
		double Iin = TransferInv(beta * mSeqArray[i], P2);
		vrestMod[i] = a * Iin * Iin + b * Iin + c;
	}


	for (int i = 0; i < NoLNs; i++) {
		theLNs[i]->SetVrest(vrestMod[i]);
		}

	//cout << endl << "All the sme rest p[oetials" << endl;

	for (int j = 0; j < NoLNs; j++) {
	for (int i = 0; i < NoLNs; i++) {


			if (i != j)
				theSynapses[j * NoLNs + i]->set_gsyn(mWeights[i][j]);

		}

	}


	vrestArray = vrestMod;
}


void neuroSynNetworkRate::ScaleWeights(double factor) {


	for (int j = 0; j < NoLNs; j++) {
	for (int i = 0; i < NoLNs; i++) {


			if (i != j)
				theSynapses[j * NoLNs + i]->set_gsyn(mWeights[i][j]*factor);

		}

	}
}



void neuroSynNetworkRate::PrintWeights()
{
	cout << "The  weighst are: (rate)" << endl;
	for (int j = 0; j < NoLNs; j++) {
		for (int i = 0; i < NoLNs; i++) {
			cout << mWeights[i][j] * (vrestMod[j] - PARAMS_LN[12]) << " ";
		}
		cout << endl;
	}
	cout << endl;


}


Array1D<double> neuroSynNetworkRate::SetCC() {

	Array1D<double> ccBias(NoLNs, 0.0);

	double beta = PARAMS_LN[15];

	double Parameters[NEUROSYN_PNO];

	for (int i = 0; i < NEUROSYN_PNO; i++)
		Parameters[i] = PARAMS_LN[i];

	vector<double> summer(NoLNs, 0);

	for (int j = 0; j < NoLNs; j++)
		for (int i = 0; i < NoLNs; i++)
			summer[j] = summer[j] + mWeights[i][j] * (vrestMod[j]
					- PARAMS_LN[12]) * mSeqArray[i];

	for (int i = 0; i < NoLNs; i++) {
		Parameters[11] = TransferInv(beta * mSeqArray[i], Parameters)
				- summer[i];
		ccBias[i] = Parameters[11];
		theLNs[i]->set_p(Parameters);
	}

	mCCBias = ccBias;
	return ccBias;
}

void neuroSynNetworkRate::PrintThetaValues()
{
	cout << "The theta values are (rate):" << endl;
	for (int i = 0; i < NoLNs; i++) {


			cout <<mCCBias[i]<< " ";

		}
		cout << endl;
}


Array2D<double> neuroSynNetworkRate::SetCritical(double epsilon) {

	Array2D<double> modWeigths(NoLNs, NoLNs, 0.0);

	for (int i = 0; i < NoLNs; i++) {
		for (int j = 0; j < NoLNs; j++) {
			modWeigths[i][j] = mWeights[i][j] * (vrestMod[j] - PARAMS_LN[12]);
		}
	}

	double maxeig = GetMaxEig();

	for (int i = 0; i < NoLNs; i++)
		for (int j = 0; j < NoLNs; j++)
			mWeights[i][j] = mWeights[i][j] * PARAMS_LN[15] / (maxeig)
					* (epsilon);

	SetWeights(mWeights, mSeqArray);
	return mWeights;
}

double neuroSynNetworkRate::GetMaxEig() {
	Array2D<double> lWeigths(NoLNs, NoLNs, 0.0);
	for (int i = 0; i < NoLNs; i++) {
		for (int j = 0; j < NoLNs; j++) {
			double gammaDeriv = TransferDeriv(TransferInv(PARAMS_LN[15]
					* mSeqArray[j], PARAMS_LN), PARAMS_LN);
			lWeigths[i][j] = mWeights[i][j] * (vrestMod[j] - PARAMS_LN[12])
					* gammaDeriv;

		//	cout  << gammaDeriv<< endl;
		}
	}

	JAMA::Eigenvalue<double> eigenDecomposition(lWeigths);
	TNT::Array1D<double> eigenValues;
	eigenDecomposition.getRealEigenvalues(eigenValues);

	double maxeig = -100000.0;
	for (int i = 0; i < NoLNs; i++) {
		if (eigenValues[i] > maxeig)
			maxeig = eigenValues[i];
	}

	return maxeig;
}

void neuroSynNetworkRate::PrintMaxEig()
{
	double lMaxEig = GetMaxEig();

	cout << " The max eig is:" << lMaxEig << endl;

}

void neuroSynNetworkRate::init() {

	dt = 0.0001;
	dtx = 0.0;
	N = 0.0;

	int counter = 0;
	int counter2 = 0;
	for (niter = neurs.begin(); niter != neurs.end(); niter++) {
		double temp[1];

      if (counter  < NoLNs) {
			temp[0] = mSeqArray[counter2];
			counter2++;
		}
		counter = counter + 1;
	}

	for (int i = 0; i < NoDirectLNInput; i++)
		directLNInput[i]->set_I(0.0);
}

Array2D<double> neuroSynNetworkRate::getJacobian(void) {

	Array2D<double> jacobian(NoLNs, NoLNs, 0.0);
	for (int i = 0; i < NoLNs; i++) {
		for (int j = 0; j < NoLNs; j++) {

			double gammaDeriv = TransferDeriv(TransferInv(PARAMS_LN[15]
					* mSeqArray[j], PARAMS_LN), PARAMS_LN);
			jacobian[i][j] = mWeights[i][j] * (vrestMod[j] - PARAMS_LN[12])
					* gammaDeriv;
			if (i == j) {
				jacobian[i][j] = -PARAMS_LN[15];
			}
		}
	}

	JAMA::Eigenvalue<double> eigenDecomposition(jacobian);
	TNT::Array1D<double> eigenValues;
	eigenDecomposition.getRealEigenvalues(eigenValues);

	cout << "Its eignevalue is:";
	double marker = -100000.0;
	for (int i = 0; i < NoLNs; i++) {
		if (eigenValues[i] > marker)
			marker = eigenValues[i];
	}
	cout << marker << " ";
	cout << endl;
	cout << flush;

	return jacobian;
}

Array1D<double> neuroSynNetworkRate::run(double tme, double inputCurr,Array1D<int> inputStart) {





	Array1D<double> endPoints(NoLNs * 2, 0.0);

	stringstream name;
	char thename[80];
	ofstream NSRateDataN, NSRateDataS;

	NSRateDataN.precision(10);
	NSRateDataS.precision(10);
	name.clear();

	int fileInt = int(inputCurr);
	if(plotInc)
	name << "NS"<< 100*percentCritical  << "RateDataF" << int(10000000* inputCurr ) << ".dat" << ends;
	else
	name << "RateF"<< "N"  << "fs"  <<  ".dat" << ends;

	//	name << "NS"  << "RateDataN"   << ".dat" << ends;
	name >> thename;

	if (OutDat)
		NSRateDataN.open(thename);

	name.clear();
	if(plotInc)
	name << "NS"<< 100*percentCritical  << "RateDataS" << int(10000000* inputCurr ) << ".dat" << ends;
	else
	name << "RateS"<< "N"  << "fs"  <<  ".dat" << ends;
	//	name << "NS"   << "RateDataS" << ".dat" << ends;
	name >> thename;
	if (OutDat)
		NSRateDataS.open(thename);

	double *tmp;

	x[0] = 0;

	bool once = true;
	while (x[0] < tme) {

		for (int i = 0; i < NoLNs; i++) {
			//		cerr << theLNs[i]->E(x) << " ";
			if (isnan(theLNs[i]->E(x))) {
				//		cerr << "nan encountered!" << endl;
				exit(1);
			}
		}
		//		cerr << endl;

		double tdt = 0.0;
		while (abs(tdt - 0.1) > 1e-9) {


				for (int i = 0; i < NoDirectLNInput; i++){
					if ((x[0] > IMPULSESTART+inputStart[i]) && (x[0] < (IMPULSESTART+inputStart[i] + 1.0)))
					directLNInput[i]->set_I(inputCurr);
			}

				for (int i = 0; i < NoDirectLNInput; i++){

			if (x[0] > IMPULSESTART+inputStart[i] + IMPULSEDUR && (x[0] < (IMPULSESTART+inputStart[i]
					+ IMPULSEDUR + 1.0)))
					directLNInput[i]->set_I(0.0);
			}

			dt = min(dt, 0.1 - tdt);
			dtx = machine->integrate(x, xn, model, dt);

			dtx = min(dtx, 2.0 * dt);
			tmp = x;
			x = xn;
			xn = tmp;
			tdt += dt;
			dt = dtx;
		}

		//	cout << x[0] << endl;

		if (OutDat) {
			NSRateDataN << x[0];


			for (int i = 0; i < NoLNs; i++) {
				NSRateDataN << " " << theLNs[i]->F(x);
			}



			NSRateDataN << endl;

			NSRateDataS << x[0];



			for (int i = 0; i < NoLNs; i++) {
				NSRateDataS << " " << theLNs[i]->S(x);
			}


			NSRateDataS << endl;
		}

		if (once) {
			for (int i = 0; i < NoLNs; i++) {
				endPoints[i] = theLNs[i]->F(x);
			}
			once = false;
		}

	}

	if (OutDat) {
		NSRateDataN.close();
		NSRateDataS.close();
	}

	for (int i = NoLNs; i < 2* NoLNs ; i++) {
		endPoints[i] = theLNs[i - NoLNs]->F(x);
	}

	return endPoints;

}

double neuroSynNetworkRate::TransferInv(double in, double *params) {

	double alpha = params[14];
	double beta = params[15];
	double tr = params[16];

	double lF;
	lF = in / alpha / tr;

	double out;

	if (theTrout)
		out = (pow((2* lF + A2 * B), 2) - pow(A2 * B, 2)) / 4 / A2 + Ic;
	else
		out = (lF - C) / M;

	out = max(out, 0.0);

	return out;

}


double neuroSynNetworkRate::TransferDeriv(double in, double *params) {

	return (Transfer(in + 0.00005, PARAMS_LN) - Transfer(in - 0.00005,
			PARAMS_LN)) / 0.0001;

}



double neuroSynNetworkRate::Transfer(double in, double *params) {

	double alpha = params[14];
	double beta = params[15];
	double tr = params[16];

	double currentIn = max(in, 0.0);

	double lF;

	if (theTrout) {
		lF = (sqrt(4 * A2 * (currentIn) + pow((A2 * B), 2)) - A2 * B) / 2;
	} else
		lF = M * currentIn + C;

	lF = max(lF, 0.0);

	double synConst = PARAMS_LN[14]
			* (exp(PARAMS_LN[15] * PARAMS_LN[16]) - 1.0);

	return alpha * tr * lF;

}


#endif

