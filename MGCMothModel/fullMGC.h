/*--------------------------------------------------------------------------
 Author: Christopher L Buckley

 Institute: Centre for Computational Neuroscience and Robotics
 University of Sussex
 FcondNetworkmer, Brighton BN1 9QJ, UK

 email to:  c.l.buckley@sussex.ac.uk

 version: 2009-11-02
 --------------------------------------------------------------------------*/

#ifndef fullMGC_H
#define fullMGC_H

class fullMGC {
public:
	list<neuron *> neurs;
	list<synapse *> syns;
	list<neuron *>::iterator niter;
	list<synapse *>::iterator siter;

	NeuroSynAdapt **theNeurons;
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
	fullMGC();
	~fullMGC();

	void generateNetwork();
	void generateConnect();
	void init();
	void enable();
	Array1D<double> run(double tme, double inputCurr, Array1D<int> inputStart);
	void SetWeights(Array2D<double> weights, Array1D<double> SeqArray);
	void SetCC(Array1D<double> biasCC);
	void ScaleWeights(double factor);
	void ScaleBias(double factor);

};

fullMGC::fullMGC() {
	dt = 0.0001;
	dtx = 0.0;
	N = 0;

	directInput = new DCInput*[NoORNs];
	theNeurons = new NeuroSynAdapt*[NoTotal];
	theORNs = new PoissonRateNeuron*[NoORNs];

	//Set up LNs and PNs
	for (int i = 0; i < NoTotal; i++) {
		if (i < NoFast)
					theNeurons[i] = new NeuroSynAdapt(i, PARAMS_LNFAST);

				if(i>=NoFast && i<NoLNs)
					theNeurons[i] = new NeuroSynAdapt(i, PARAMS_LN);

		        if(i>=NoLNs)
					theNeurons[i] = new NeuroSynAdapt(i, PARAMS_PN);

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

fullMGC::~fullMGC() {
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

void fullMGC::enable() {
	model = new NeuronModel(&neurs, &syns, N, cerr);
	x = new double[N];
	xn = new double[N];
	machine = new rk65n(N, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);
}

void fullMGC::ScaleWeights(double factor) {

	for (int j = 0; j < NoTotal; j++) {
		for (int i = 0; i < NoTotal; i++) {
			if (i != j)
				theSynapses[j * NoTotal + i]->set_gsyn(mWeights[i][j] * factor);
		}
	}
}

void fullMGC::ScaleBias(double factor) {

	for (int i = 0; i < NoTotal; i++) {
		theNeurons[i]->p[11] = (theNeurons[i]->p[11] * factor);
	}
}

void fullMGC::SetWeights(Array2D<double> weights, Array1D<double> SeqArray) {

	mWeights = weights;
	mSeqArray = SeqArray;

	for (int i = 0; i < NoTotal; i++)
		theNeurons[i]->p[18] = SeqArray[i];

	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {

			if (i != j) {
				theSynapses[j * NoTotal + i]->set_gsyn(mWeights[i][j]);
			}
		}
	}

}

void fullMGC::SetCC(Array1D<double> biasCC) {
	for (int i = 0; i < NoTotal; i++)
		theNeurons[i]->p[11] = biasCC[i];
}

void fullMGC::init() {
	dt = 0.0001;
	dtx = 0.0;
	N = 0;

	double Initvalue[NEUROSYNADAPT_IVARNO];
	for (niter = neurs.begin(); niter != neurs.end(); niter++) {

		for (int i = 0; i < NEUROSYNADAPT_IVARNO; i++)
			Initvalue[i] = NEUROSYNADAPT_INIVARS[i] * RG.n() - 0.5 * RG.n();

		Initvalue[6] = 0.0;
		Initvalue[7] =0.0;
		Initvalue[5] = RG.n();
		(*niter)->init(x, Initvalue);

	}



	for (int i = 0; i < NoORNs; i++)
		directInput[i]->set_I(0.0);

}

Array1D<double> fullMGC::run(double tme, double inputCurr,
		Array1D<int> inputStart) {



	Array1D<double> endPoints(2 * NoTotal, 0.0);
	vector<double> spike_history;
	stringstream name;
	char thename[80];
	ofstream NSDataN, NSDataS, NSDataIysn, NSDataM, NSDataCa, NSDataTheta,
			NSDataOM, NSDataOS, NSDataON;

	NSDataN.precision(10);
	NSDataS.precision(10);
	NSDataIysn.precision(10);
	NSDataM.precision(10);
	NSDataCa.precision(10);
	NSDataTheta.precision(10);
	NSDataOM.precision(10);
	NSDataOS.precision(10);
	NSDataON.precision(10);
	name.clear();

	double *tmp;

	x[0] = 0;

	double factor;

	if (doneFileCreate) {
		if (Manual) {
		/*	name << globalName << "CondF.dat" << ends;
			name >> thename;
			NSDataN.open(thename);
			name.clear();
*/
			name << globalName << "CondM.dat" << ends;
			name >> thename;
			NSDataM.open(thename);
			name.clear();
/*
			name << globalName << "CondS.dat" << ends;
			name >> thename;
			NSDataS.open(thename);
			name.clear();

			name << globalName << "CondCa.dat" << ends;
			name >> thename;
			NSDataCa.open(thename);
			name.clear();

			name << globalName << "CondTheta.dat" << ends;
			name >> thename;
			NSDataTheta.open(thename);
			name.clear();

			name << globalName << "CondOS.dat" << ends;
			name >> thename;
			NSDataOS.open(thename);
			name.clear();

			name << globalName << "CondON.dat" << ends;
			name >> thename;
			NSDataON.open(thename);
			name.clear();

			name << globalName << "CondOM.dat" << ends;
			name >> thename;
			NSDataOM.open(thename);
			name.clear();
			*/
		} else {
	/*		name << "CondF.dat" << ends;
			name >> thename;
			NSDataN.open(thename);
			name.clear();
*/
			name << "CondM.dat" << ends;
			name >> thename;
			NSDataM.open(thename);
			name.clear();
/*
			name << "CondS.dat" << ends;
			name >> thename;
			NSDataS.open(thename);
			name.clear();

			name << "CondCa.dat" << ends;
			name >> thename;
			NSDataCa.open(thename);
			name.clear();

			name << "CondTheta.dat" << ends;
			name >> thename;
			NSDataTheta.open(thename);
			name.clear();

			name << "CondOS.dat" << ends;
			name >> thename;
			NSDataOS.open(thename);
			name.clear();

			name << "CondON.dat" << ends;
			name >> thename;
			NSDataON.open(thename);
			name.clear();

			name << "CondOM.dat" << ends;
			name >> thename;
			NSDataOM.open(thename);
			name.clear();
			*/
		}

	}

	while (x[0] < tme) {

		if (Gradual) {

			if (int(x[0] * 10) % int(scaleBin * 10) == 0 && x[0] < scaleTime) {

				factor = (scaleBin + x[0]) / scaleTime;
				ScaleWeights(0.5 + 0.3 * factor);
				//		cout << (0.5 + 0.4 * factor) * percentCritical << endl;

			}

			if (x[0] > scaleTime && x[0] < TURNOFFADAPT)
				CaAdapt = true;
			else
				CaAdapt = false;

			if (int((x[0] - scaleStart2) * 10) % int(scaleBin * 10) == 0
					&& x[0] > scaleStart2 && x[0] < (scaleStart2 + scaleTime2)) {

				factor = (scaleBin + x[0] - scaleStart2) / scaleTime2;
				ScaleWeights(0.8 + 0.2 * factor);
				//		cout << (0.9  + 0.1*factor) * percentCritical << endl;

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

		double tdt = 0.0;

		for (int i = 0; i < NoTotal; i++) {
			if (isnan(theNeurons[i]->E(x))) {
				exit(1);
			}

			if (theNeurons[i]->S(x) < 0.001 && x[0]
					> (scaleStart2 + scaleTime2) && x[0] < IMPULSESTART
					&& CreateFile) {
				//	cout << "I touched the bottom so I am exiting" <<endl;
				//	exit(1);
			}
		}

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

		for (int i = 0; i < NoTotal; i++)
			theNeurons[i]->spike_detect(x);

	/*	for (int i = 0; i < NoORNs; i++) {
				theORNs[i]->validate_E(x, tdt);
				theORNs[i]->step();
			}
*/

		//out to screen progress
		if (int(x[0] * 10) % 5000 == 0)
				cout  << x[0] <<endl;


	//	for (int i = 0; i < 18; i++)
	//			cout << x[i] << " ";

	//	cout << endl;
		if (int(x[0] * 10) % 10 == 0)
				{
		if (doneFileCreate && x[0] > OutTime) {
		//	NSDataN << x[0];
		//	NSDataCa << x[0];
		//	NSDataTheta << x[0];
	//		NSDataS << x[0];
			NSDataM << x[0];
		//	NSDataON << x[0];
		//	NSDataOS << x[0];
		//	NSDataOM << x[0];

			for (int i = 0; i < NoORNs; i++) {

		//		NSDataON << " " << theORNs[i]->F(x);
		//		NSDataOS << " " << theORNs[i]->S(x);
		//		NSDataOM << " " << 1;

			}

			for (int i = 0; i < NoTotal; i++) {
				double spiker = 0;
				if (theNeurons[i]->start_spiking)
					spiker = 1;
		//		NSDataN << " " << spiker;
			//	NSDataCa << " " << theNeurons[i]->Ca(x);
			//	NSDataTheta << " " << theNeurons[i]->Theta(x);
		//		NSDataS << " " << theNeurons[i]->S(x);
				NSDataM << " " << theNeurons[i]->E(x);
			}

	//		NSDataN << endl;
		//	NSDataCa << endl;
		//	NSDataTheta << endl;
	//		NSDataS << endl;
			NSDataM << endl;
		//	NSDataON << endl;
		//	NSDataOS << endl;
		//	NSDataOM << endl;
		}
				}

		if (x[0] > (IMPULSESTART - 800) && x[0] < IMPULSESTART) {

			for (int i = 0; i < NoTotal; i++) {
				double spiker = 0;

				if (theNeurons[i]->start_spiking)
					spiker = 1.0;

				endPoints[i] = endPoints[i] + spiker;
			}

		}

		if (x[0] > (THETIME - 800)) {

			for (int i = NoTotal; i < 2 * NoTotal; i++) {
				double spiker = 0;

				if (theNeurons[i - NoTotal]->start_spiking)
					spiker = 1.0;

				endPoints[i] = endPoints[i] + spiker;
			}

		}

	}

	if (doneFileCreate) {
		/*NSDataON.close();
		NSDataOS.close();
		NSDataM.close();
		NSDataN.close();
		NSDataS.close();
		NSDataCa.close();
		NSDataTheta.close();
		*/
		NSDataM.close();
	}

	for (int i = 0; i < 2 * NoTotal; i++)
		endPoints[i] = endPoints[i] / 800.0;

	return endPoints;

}

#endif
