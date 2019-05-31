/*--------------------------------------------------------------------------
 Author: Christopher L Buckley

 Institute: Centre for Computational Neuroscience and Robotics
 University of Sussex
 FcondNetworkmer, Brighton BN1 9QJ, UK

 email to:  c.l.buckley@sussex.ac.uk

 version: 2009-11-02
 --------------------------------------------------------------------------*/

#ifndef fullMGCRate_H
#define fullMGCRate_H


//w[i][j] source i  and target j

class fullMGCRate {
public:
	list<neuron *> neurs;
	list<synapse *> syns;
	list<neuron *>::iterator niter;
	list<synapse *>::iterator siter;

	NeuroSynRate **theLNs;
	NeuroSynRate **theORNs;
	Emptysynapse **theSynapses;

	DCInput **directLNInput;
	DCInput **directORNInput;
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

	fullMGCRate();
	~fullMGCRate();

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
	Array2D<double> SetCritical(double epsilon);
	Array1D<double> SetCC();
	double GetMaxEig();
};

fullMGCRate::fullMGCRate() {
	dt = 0.0001;
	dtx = 0.0;
	N = 0.0;

	directLNInput = new DCInput*[NoTotal];
	directORNInput = new DCInput*[NoORNs];
	theLNs = new NeuroSynRate*[NoTotal];
	theORNs = new NeuroSynRate*[NoORNs];


	//Set up network memory and space


	for (int i = 0; i < NoTotal; i++) {
		theLNs[i] = new NeuroSynRate(i, PARAMS_LN);
		neurs.push_back(theLNs[i]);
		directLNInput[i] = new DCInput(theLNs[i], 1);
	}

	//set up connections
	theSynapses = new Emptysynapse*[NoTotal * NoTotal];
	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {

			if (i != j) {
				theSynapses[j * NoTotal + i] = new Emptysynapse(theLNs[i],
						theLNs[j], 0.0);
			}
		}
	}


	for (int i = 0; i < NoORNs; i++) {
			theORNs[i] = new NeuroSynRate(i, PARAMS_ORN);
			neurs.push_back(theORNs[i]);
			directORNInput[i] = new DCInput(theORNs[i], 0.0);
		}

	for (int i = 0; i < NoORNs; i++) {
			Emptysynapse *aSynapse;
			aSynapse = new Emptysynapse(theORNs[i], theLNs[i], -ORN2LN);

		}

	//enable integaror
	enable();
}

fullMGCRate::~fullMGCRate() {
	list<neuron *>::iterator i;
	list<synapse *>::iterator j;
	for (i = neurs.begin(); i != neurs.end(); i++) {
		delete *i;
	}
	for (j = syns.begin(); j != syns.end(); j++) {
		delete *j;
	}

	delete[] theLNs;
	delete[] theORNs;
	delete[] directLNInput;
	delete[] directORNInput;


}

void fullMGCRate::enable() {
	model = new NeuronModel(&neurs, &syns, N, cerr);
	x = new double[N];
	xn = new double[N];
	machine = new rk65n(N, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);
}

void fullMGCRate::SetWeights(Array2D<double> weights,
		Array1D<double> SeqArray) {

	mWeights = weights;
	mSeqArray = SeqArray;

	double beta = PARAMS_LN[15];


	//this code alters the wigtsh absed on the rest potetial.

	Array1D<double> setter(NoTotal, 0.0);

	vrestMod = setter;

	double P2[NEUROSYN_PNO];

	for (int i = 0; i < NEUROSYN_PNO; i++)
		P2[i] = PARAMS_LN[i];

	for (int i = 0; i < NoTotal; i++) {
		double Iin = TransferInv(beta * mSeqArray[i], P2);
		vrestMod[i] = a * Iin * Iin + b * Iin + c;
	}


	for (int i = 0; i < NoTotal; i++) {
		theLNs[i]->SetVrest(vrestMod[i]);
		}

	//cout << endl << "All the sme rest p[oetials" << endl;

	for (int j = 0; j < NoTotal; j++) {
	for (int i = 0; i < NoTotal; i++) {


			if (i != j)
				theSynapses[j * NoTotal + i]->set_gsyn(mWeights[i][j]);

		}

	}


}

void fullMGCRate::PrintWeights()
{
	cout << "The  weighst are: (rate)" << endl;
	for (int j = 0; j < NoTotal; j++) {
		for (int i = 0; i < NoTotal; i++) {
			cout << mWeights[i][j] * (vrestMod[j] - PARAMS_LN[12]) << " ";
		}
		cout << endl;
	}
	cout << endl;


}


Array1D<double> fullMGCRate::SetCC() {

	Array1D<double> ccBias(NoTotal, 0.0);

	double beta = PARAMS_LN[15];

	double Parameters[NEUROSYN_PNO];

	for (int i = 0; i < NEUROSYN_PNO; i++)
		Parameters[i] = PARAMS_LN[i];

	vector<double> summer(NoTotal, 0);

	for (int j = 0; j < NoTotal; j++)
		for (int i = 0; i < NoTotal; i++)
			summer[j] = summer[j] + mWeights[i][j] * (vrestMod[j]
					- PARAMS_LN[12]) * mSeqArray[i];

	for (int i = 0; i < NoTotal; i++) {
		Parameters[11] = TransferInv(beta * mSeqArray[i], Parameters)
				- summer[i];
		ccBias[i] = Parameters[11];
		theLNs[i]->set_p(Parameters);
	}

	mCCBias = ccBias;
	return ccBias;
}

void fullMGCRate::PrintThetaValues()
{
	cout << "The theta values are (rate):" << endl;
	for (int i = 0; i < NoTotal; i++) {


			cout <<mCCBias[i]<< " ";

		}
		cout << endl;
}


Array2D<double> fullMGCRate::SetCritical(double epsilon) {

	Array2D<double> modWeigths(NoTotal, NoTotal, 0.0);

	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {
			modWeigths[i][j] = mWeights[i][j] * (vrestMod[j] - PARAMS_LN[12]);
		}
	}

	double maxeig = GetMaxEig();

	for (int i = 0; i < NoTotal; i++)
		for (int j = 0; j < NoTotal; j++)
			mWeights[i][j] = mWeights[i][j] * PARAMS_LN[15] / (maxeig)
					* (epsilon);

	SetWeights(mWeights, mSeqArray);
	return mWeights;
}

double fullMGCRate::GetMaxEig() {
	Array2D<double> lWeigths(NoTotal, NoTotal, 0.0);
	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {
			double gammaDeriv = TransferDeriv(TransferInv(PARAMS_LN[15]
					* mSeqArray[j], PARAMS_LN), PARAMS_LN);
			lWeigths[i][j] = mWeights[i][j] * (vrestMod[j] - PARAMS_LN[12])
					* gammaDeriv;
		}
	}

	JAMA::Eigenvalue<double> eigenDecomposition(lWeigths);
	TNT::Array1D<double> eigenValues;
	eigenDecomposition.getRealEigenvalues(eigenValues);

	double maxeig = -100000.0;
	for (int i = 0; i < NoTotal; i++) {
		if (eigenValues[i] > maxeig)
			maxeig = eigenValues[i];
	}

	return maxeig;
}

void fullMGCRate::PrintMaxEig()
{
	double lMaxEig = GetMaxEig();

	cout << " The max eig is:" << lMaxEig << endl;

}

void fullMGCRate::init() {

	dt = 0.0001;
	dtx = 0.0;
	N = 0.0;

	int counter = 0;
	int counter2 = 0;
	for (niter = neurs.begin(); niter != neurs.end(); niter++) {
		double temp[1];


		 if (counter < NoTotal) {
			temp[0] = mSeqArray[counter2];
			counter2++;
		}
		 else
		 {
			 temp[0] = 0;
		 }

		(*niter)->init(x, temp);

		counter = counter + 1;
	}

	for (int i = 0; i < NoTotal; i++)
		directLNInput[i]->set_I(0.0);
	for (int i = 0; i < NoORNs; i++)
		directORNInput[i]->set_I(0.0);
}

Array2D<double> fullMGCRate::getJacobian(void) {

	Array2D<double> jacobian(NoTotal, NoTotal, 0.0);
	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {

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

//	cout << "Its eignevalue is:";
	double marker = -100000.0;
	for (int i = 0; i < NoTotal; i++) {
		if (eigenValues[i] > marker)
			marker = eigenValues[i];
	}
//	cout << marker << " ";
//	cout << endl;
//	cout << flush;

	return jacobian;
}

Array1D<double> fullMGCRate::run(double tme, double inputCurr,Array1D<int> inputStart) {





	Array1D<double> endPoints(NoTotal * 2, 0.0);

	stringstream name;
	char thename[80];
	ofstream NSRateDataN, NSRateDataS;

	NSRateDataN.precision(10);
	NSRateDataS.precision(10);
	name.clear();


	if(plotInc)
	name << globalName << "RateDataF" << int(10000000* inputCurr ) << ".dat" << ends;
	else
	{
		if(Manual)
			name << globalName << "RateF"<< "N"  << "fs"  <<  ".dat" << ends;
		else
			name << "RateF"<< "N"  << "fs"  <<  ".dat" << ends;
	}

	//	name << "NS"  << "RateDataN"   << ".dat" << ends;
	name >> thename;

	if (OutDat)
		NSRateDataN.open(thename);

	name.clear();
	if(plotInc)
	name << globalName  << "RateDataS" << int(10000000* inputCurr ) << ".dat" << ends;
	else
	{
		if(Manual)
			name << globalName << "RateS"<< "N"  << "fs"  <<  ".dat" << ends;
		else
			name << "RateS"<< "N"  << "fs"  <<  ".dat" << ends;
	}
	//	name << "NS"   << "RateDataS" << ".dat" << ends;
	name >> thename;
	if (OutDat)
		NSRateDataS.open(thename);

	double *tmp;

	x[0] = 0;

	bool once = true;
	while (x[0] < tme) {

		for (int i = 0; i < NoTotal; i++) {
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
					if ((x[0] > IMPULSESTART+inputStart[i]) && (x[0] < (IMPULSESTART+inputStart[i] + 1.0))){
					directORNInput[i]->set_I(inputCurr);
					}
			}

				for (int i = 0; i < NoDirectLNInput; i++){

			if (x[0] > IMPULSESTART+inputStart[i] + IMPULSEDUR && (x[0] < (IMPULSESTART+inputStart[i]
					+ IMPULSEDUR + 1.0)))
				directORNInput[i]->set_I(0.0);
			}



				for (int i = NoLNs; i < NoTotal; i++){
					if ((x[0] > IMPULSESTART) && (x[0] < (IMPULSESTART + 1.0)))
						directORNInput[i]->set_I(inputCurr*ORN2PN);
			}

				for (int i = NoLNs; i < NoTotal; i++){

			if (x[0] > IMPULSESTART + IMPULSEDUR && (x[0] < (IMPULSESTART
					+ IMPULSEDUR + 1.0)))
				directORNInput[i]->set_I(0.0);
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



		if (OutDat) {
			NSRateDataN << x[0];

			for (int i = 0; i < NoTotal; i++) {
				NSRateDataN << " " << theLNs[i]->F(x);
			}

			for (int i = 0; i < NoORNs; i++)
							NSRateDataN << " " << theORNs[i]->F(x);

			NSRateDataN << endl;

			NSRateDataS << x[0];

			for (int i = 0; i < NoTotal; i++) {
				NSRateDataS << " " << theLNs[i]->S(x);
			}

			for (int i = 0; i < NoORNs; i++)
				NSRateDataS << " " << theORNs[i]->S(x);



			NSRateDataS << endl;
		}

		if (once) {
			for (int i = 0; i < NoTotal; i++) {
				endPoints[i] = theLNs[i]->F(x);
			}
			once = false;
		}

	}

	if (OutDat) {
		NSRateDataN.close();
		NSRateDataS.close();
	}

	for (int i = NoTotal; i < 2* NoTotal ; i++) {
		endPoints[i] = theLNs[i - NoTotal]->F(x);
	}

	return endPoints;

}


double fullMGCRate::TransferInv(double in, double *params) {

	double alpha = params[14];
	double tr = params[16];


	double lF;
	lF =  in/alpha/tr;

	double out;


	if (theTrout)
		out = (pow((2* lF + A2 * B), 2) - pow(A2 * B, 2)) / 4 / A2 +Ic;
	else
		out = (lF - C) / M;


	 out = max(out,0.0);



	return out;

}


double fullMGCRate::TransferDeriv(double in, double *params) {

	return (Transfer(in + 0.00005, PARAMS_LN) - Transfer(in - 0.00005,
			PARAMS_LN)) / 0.0001;

}



double fullMGCRate::Transfer(double in, double *params) {

	double alpha = params[14];
	double tr = params[16];

	double currentIn = max(in, 0.0);

	double lF;

	if(theTrout)
		{
		 lF = (sqrt(4 * A2 * (currentIn-Ic) + pow((A2*B), 2)) - A2
				*B) / 2;
		}
		else
			lF = M * currentIn + C;

	lF = max(lF, 0.0);
	return alpha*tr*lF;
}


#endif

