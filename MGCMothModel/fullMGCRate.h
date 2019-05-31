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

	NeuroSynRate **theNeurons;
	Emptysynapse **theSynapses;

	PoissonRateNeuron **theORNs;
	Emptysynapse **theORNSynapses;

	DCInput **directInput;
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

	Array1D<double> run(double tme, double inputCurr, Array1D<int> inputStart);
	double TransferInv(double in, double *params);
	double Transfer(double in, double *params);
	double TransferDeriv(double in, double *params);
	Array2D<double> getJacobian(void);
	void SetWeights(Array2D<double> weights, Array1D<double> SeqArray);
	Array2D<double> SetCritical(double epsilon);
	Array1D<double> SetCC();
	void ScaleWeights(double factor);
	void ScaleBias(double factor);
	double GetMaxEig();
};

fullMGCRate::fullMGCRate() {
	dt = 0.0001;
	dtx = 0.0;
	N = 0;

	directInput = new DCInput*[NoORNs];
	theNeurons = new NeuroSynRate*[NoTotal];
	theORNs = new PoissonRateNeuron*[NoORNs];

	//Set up LNs and PNs
	for (int i = 0; i < NoTotal; i++) {
		if (i < NoFast)
					theNeurons[i] = new NeuroSynRate(i, PARAMS_LNFAST);

				if(i>=NoFast && i<NoLNs)
					theNeurons[i] = new NeuroSynRate(i, PARAMS_LN);

		        if(i>=NoLNs)
					theNeurons[i] = new NeuroSynRate(i, PARAMS_PN);

		neurs.push_back(theNeurons[i]);
	//	directInput[i] = new DCInput(theNeurons[i], 1.0);
	}

	//The ORNs
	for (int i = 0; i < NoORNs; i++) {
		theORNs[i] = new PoissonRateNeuron(i, PARAMS_ORN);
		neurs.push_back(theORNs[i]);
		directInput[i] = new DCInput(theORNs[i], 1.0);
	}

	//set up  AL connections
	theSynapses = new Emptysynapse*[NoTotal * NoTotal];
	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {
			if (i != j) {
				theSynapses[j * NoTotal + i] = new Emptysynapse(theNeurons[i],
						theNeurons[j], 0.0);
			}
		}
	}

	//set up  ORN-AL connections
	theORNSynapses = new Emptysynapse*[NoORNs + NoORNs * NoPNs];
	for (int i = 0; i < NoORNs; i++) {
		theORNSynapses[i] = new Emptysynapse(theORNs[i], theNeurons[i+stimNum],
				-ORN2LNWeight);
	}

	for (int i = 0; i < NoORNs; i++) {
		for (int j = NoLNs; j < NoTotal; j++) {
			theORNSynapses[NoORNs + (j - NoLNs) * NoORNs + i]
					= new Emptysynapse(theORNs[i], theNeurons[j],
							-ORN2PNWeight);
		}
	}

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

	delete[] theNeurons;
	delete[] theORNs;
	delete[] directInput;

}

void fullMGCRate::enable() {
	model = new NeuronModel(&neurs, &syns, N, cerr);
	x = new double[N];
	xn = new double[N];
	machine = new rk65n(N, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);
}
void fullMGCRate::ScaleWeights(double factor) {

	for (int j = 0; j < NoTotal; j++) {
		for (int i = 0; i < NoTotal; i++) {

			if (i != j)
				theSynapses[j * NoTotal + i]->set_gsyn(mWeights[i][j] * factor);

		}

	}
}

void fullMGCRate::ScaleBias(double factor) {
	for (int i = 0; i < NoTotal; i++)
		theNeurons[i]->p[11] = (theNeurons[i]->p[11] * factor);
}

void fullMGCRate::SetWeights(Array2D<double> weights, Array1D<double> SeqArray) {

	mWeights = weights;
	mSeqArray = SeqArray;

	//this code alters the wigtsh absed on the rest potetial.

	Array1D<double> setter(NoTotal, 0.0);

	vrestMod = setter;
	for (int i = 0; i < NoTotal; i++) {
		double Iin = TransferInv(theNeurons[i]->p[15] * mSeqArray[i],
				theNeurons[i]->p);
		vrestMod[i] = a * Iin * Iin + b * Iin + c;
	}

	for (int i = 0; i < NoTotal; i++) {
		theNeurons[i]->SetVrest(vrestMod[i]);
	}

	//cout << endl << "All the sme rest p[oetials" << endl;

	for (int j = 0; j < NoTotal; j++) {
		for (int i = 0; i < NoTotal; i++) {

			if (i != j)
				theSynapses[j * NoTotal + i]->set_gsyn(mWeights[i][j]);

		}

	}

}

void fullMGCRate::PrintWeights() {
	cout << "The  weighst are: (rate)" << endl;
	for (int j = 0; j < NoTotal; j++) {
		for (int i = 0; i < NoTotal; i++) {
			cout << mWeights[i][j] * (vrestMod[j] - theNeurons[i]->p[12])
					<< " ";
		}
		cout << endl;
	}
	cout << endl;

}

Array1D<double> fullMGCRate::SetCC() {

	Array1D<double> ccBias(NoTotal, 0.0);
	vector<double> summer(NoTotal, 0);

	for (int j = 0; j < NoTotal; j++)
		for (int i = 0; i < NoTotal; i++)
			summer[j] = summer[j] + mWeights[i][j] * (vrestMod[j]
					- theNeurons[i]->p[12]) * mSeqArray[i];

	for (int i = 0; i < NoTotal; i++) {
		ccBias[i] = TransferInv(theNeurons[i]->p[15] * mSeqArray[i],
				theNeurons[i]->p) - summer[i];
		theNeurons[i]->p[11] = ccBias[i];
		theNeurons[i]->p[18] =  mSeqArray[i];
	}

	mCCBias = ccBias;
	return ccBias;
}

void fullMGCRate::PrintThetaValues() {
	cout << "The theta values are (rate):" << endl;
	for (int i = 0; i < NoTotal; i++) {

		cout << mCCBias[i] << " ";

	}
	cout << endl;
}

Array2D<double> fullMGCRate::SetCritical(double epsilon) {

	Array2D<double> modWeigths(NoTotal, NoTotal, 0.0);

	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {
			modWeigths[i][j] = mWeights[i][j] * (vrestMod[j]
					- theNeurons[i]->p[12]);
		}
	}

	double maxeig = GetMaxEig();

	for (int i = 0; i < NoTotal; i++)
		for (int j = 0; j < NoTotal; j++)
			mWeights[i][j] = mWeights[i][j] * theNeurons[i]->p[15] / (maxeig)
					* (epsilon);

	SetWeights(mWeights, mSeqArray);
	return mWeights;
}

double fullMGCRate::GetMaxEig() {
	Array2D<double> lWeigths(NoTotal, NoTotal, 0.0);
	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {
			double gammaDeriv = TransferDeriv(TransferInv(theNeurons[j]->p[15]
					* mSeqArray[j], theNeurons[j]->p), theNeurons[j]->p);
			lWeigths[i][j] = mWeights[i][j] * (vrestMod[j]
					- theNeurons[i]->p[12]) * gammaDeriv;
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

void fullMGCRate::PrintMaxEig() {
	double lMaxEig = GetMaxEig();

	cout << " The max eig is:" << lMaxEig << endl;

}

void fullMGCRate::init() {

	dt = 0.0001;
	dtx = 0.0;
	N = 0;

	int counter = 0;
	int counter2 = 0;
	for (niter = neurs.begin(); niter != neurs.end(); niter++) {
		double temp[1];

		if (counter < NoTotal) {
			temp[0] = mSeqArray[counter2];
			counter2++;
		} else {
			temp[0] = 0;
		}

		(*niter)->init(x, temp);

		counter = counter + 1;
	}

	for (int i = 0; i < NoORNs; i++)
		directInput[i]->set_I(0.0);
}

Array2D<double> fullMGCRate::getJacobian(void) {

	Array2D<double> jacobian(NoTotal, NoTotal, 0.0);
	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {

			double gammaDeriv = TransferDeriv(TransferInv(theNeurons[j]->p[15]
					* mSeqArray[j], theNeurons[j]->p), theNeurons[j]->p);
			jacobian[i][j] = mWeights[i][j] * (vrestMod[j]
					- theNeurons[i]->p[12]) * gammaDeriv;

			if (i == j) {
				jacobian[i][j] = -theNeurons[j]->p[15];
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

Array1D<double> fullMGCRate::run(double tme, double inputCurr,
		Array1D<int> inputStart) {

	Array1D<double> endPoints(NoTotal * 2, 0.0);

	stringstream name;
	char thename[80];
	ofstream NSRateDataN, NSRateDataS, NSRateDataON, NSRateDataOS;

	NSRateDataN.precision(10);
	NSRateDataS.precision(10);
	NSRateDataON.precision(10);
	NSRateDataOS.precision(10);

	name.clear();

	double *tmp;

	x[0] = 0;
	double factor;
	bool once = true;
	if (doneFileCreate) {
		if (Manual) {
			name << globalName << "RateF.dat" << ends;
			name >> thename;
			NSRateDataN.open(thename);
			name.clear();
/*
			name << globalName << "RateS.dat" << ends;
			name >> thename;
			NSRateDataS.open(thename);
			name.clear();

			name << globalName << "RateOS.dat" << ends;
			name >> thename;
			NSRateDataOS.open(thename);
			name.clear();

			name << globalName << "RateON.dat" << ends;
			name >> thename;
			NSRateDataON.open(thename);
			name.clear();
			*/

		} else {
			name << "RateF.dat" << ends;
			name >> thename;
			NSRateDataN.open(thename);
			name.clear();
/*
			name << "RateS.dat" << ends;
			name >> thename;
			NSRateDataS.open(thename);
			name.clear();

			name << "RateOS.dat" << ends;
			name >> thename;
			NSRateDataOS.open(thename);
			name.clear();

			name << "RateON.dat" << ends;
			name >> thename;
			NSRateDataON.open(thename);
			name.clear();
			*/

		}

	}

	while (x[0] < tme) {

		if (Gradual) {

			if (int(x[0] * 10) % int(scaleBin * 10) == 0 && x[0] < scaleTime) {
				factor = (scaleBin + x[0]) / scaleTime;
				ScaleWeights(0.5 + 0.4 * factor);
			}

			if (int((x[0] - scaleStart2) * 10) % int(scaleBin * 10) == 0
					&& x[0] > scaleStart2 && x[0] < (scaleStart2 + scaleTime2)) {

				factor = (scaleBin + x[0] - scaleStart2) / scaleTime2;
				ScaleWeights(0.9 + 0.1 * factor);

			}
		}

		for (int i = 0; i < NoORNs; i++) {

			if (x[0] < IMPULSESTART || x[0] > IMPULSESTART + IMPULSEDUR)
			{
				directInput[i]->set_I(0.0);
			//	directInput[NoTotal-i-1]->set_I(0.0);
			//	for (int j = NoLNs; j < NoTotal; j++)
			//						directInput[j]->set_I(0.0);
			}
			else
			{
				directInput[i]->set_I(inputCurr);

/*
				if(doPatterned)
				{
					if (int(x[0] * 10) % 6000 >3000)
						directInput[i]->set_I(inputCurr);
					else
						directInput[i]->set_I(0);

				}
				else
					directInput[i]->set_I(inputCurr);
*/
			//	for (int j = NoLNs; j < NoTotal; j++)
			//		directInput[j]->set_I(inputCurr*ORN2PNWeight);
			}
		}

		for (int i = 0; i < NoTotal; i++) {
			if (isnan(theNeurons[i]->E(x))) {
				cerr << "nan encountered!" << endl;
				exit(1);
			}
		}

		double tdt = 0.0;
		while (abs(tdt - 0.1) > 1e-9) {

			dt = min(dt, 0.1 - tdt);
			dtx = machine->integrate(x, xn, model, dt);

			dtx = min(dtx, 2.0 * dt);
			tmp = x;
			x = xn;
			xn = tmp;
			tdt += dt;
			dt = dtx;
		}

		//out to screen progress
		if (int(x[0] * 10) % 5000 == 0)
			cout << x[0] << endl;

		if (int(x[0] * 10) % 10 == 0)
		{
		if (doneFileCreate && x[0] > OutTime) {

			NSRateDataN << x[0];
	//		NSRateDataS << x[0];
	//		NSRateDataON << x[0];
	//		NSRateDataOS << x[0];

			for (int i = 0; i < NoTotal; i++) {
				NSRateDataN << " " << theNeurons[i]->F(x);
		//		NSRateDataS << " " << theNeurons[i]->S(x);
			}
			for (int i = 0; i < NoORNs; i++) {
//				NSRateDataON << " " << theORNs[i]->F(x);
//				NSRateDataOS << " " << theORNs[i]->S(x);
			}

			NSRateDataN << endl;
		//	NSRateDataS << endl;
	//		NSRateDataON << endl;
	//		NSRateDataOS << endl;

		}
		}

		if (once && x[0] > IMPULSESTART - 0.1) {
			for (int i = 0; i < NoTotal; i++) {
				endPoints[i] = theNeurons[i]->F(x);
			}
			once = false;
		}

	}

	if (doneFileCreate) {
		NSRateDataN.close();
	//	NSRateDataS.close();
	//	NSRateDataON.close();
	//	NSRateDataOS.close();
	}

	for (int i = NoTotal; i < 2 * NoTotal; i++) {
		endPoints[i] = theNeurons[i - NoTotal]->F(x);
	}

	return endPoints;

}

double fullMGCRate::TransferInv(double in, double *params) {

	double alpha = params[14];
	double tr = params[16];

	double lF;
	lF = in / alpha / tr;

	double out;


		out = (lF - C) / M;

	out = max(out, 0.0);

	return out;

}

double fullMGCRate::TransferDeriv(double in, double *params) {

	return (Transfer(in + 0.00005, params) - Transfer(in - 0.00005, params))
			/ 0.0001;

}

double fullMGCRate::Transfer(double in, double *params) {

	double alpha = params[14];
	double tr = params[16];

	double currentIn = max(in, 0.0);

	double lF;


		lF = M * currentIn + C;

	lF = max(lF, 0.0);
	return alpha * tr * lF;
}

#endif

