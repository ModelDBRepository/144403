/*--------------------------------------------------------------------------
 Author: Christopher L Buckley

 Institute: Centre for Computational Neuroscience and Robotics
 University of Sussex
 FcondNetworkmer, Brighton BN1 9QJ, UK

 email to:  c.l.buckley@sussex.ac.uk

 version: 2009-11-02
 --------------------------------------------------------------------------*/

#ifndef neuroSynNetworkDetDET_H
#define neuroSynNetworkDetDET_H

class neuroSynNetworkDet {
public:
	list<neuron *> neurs;
	list<synapse *> syns;
	list<neuron *>::iterator niter;
	list<synapse *>::iterator siter;

	NeuroSyn **theLNs;
	NeuroSyn **thePNs;
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
	neuroSynNetworkDet();
	~neuroSynNetworkDet();
	void PrintWeights();
	void generateNetwork();
	void generateConnect();
	void init();
	void enable();
	Array1D<double> run(double tme, double inputCurr, Array1D<int> inputStart);
	void SetWeights(Array2D<double> weights, Array1D<double> SeqArray);
	void ScaleWeights(double factor);
	void ScaleBias(double factor);
	void SetCC(Array1D<double> biasCC);

};

neuroSynNetworkDet::neuroSynNetworkDet() {
	dt = 0.0001;
	dtx = 0.0;
	N = 0.0;

	directLNInput = new DCInput*[NoTotal];
	theLNs = new NeuroSyn*[NoTotal];

	//Set up network memory and space


	for (int i = 0; i < NoTotal; i++) {
		theLNs[i] = new NeuroSyn(i, PARAMS_LN);
		neurs.push_back(theLNs[i]);
		directLNInput[i] = new DCInput(theLNs[i], 0);
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

	//enable integaror
	enable();
}

neuroSynNetworkDet::~neuroSynNetworkDet() {
	list<neuron *>::iterator i;
	list<synapse *>::iterator j;
	for (i = neurs.begin(); i != neurs.end(); i++) {
		delete *i;
	}
	for (j = syns.begin(); j != syns.end(); j++) {
	}

	delete[] theLNs;
	delete[] directLNInput;
}

void neuroSynNetworkDet::PrintWeights() {
	cout << "The  weighst are: (rate)" << endl;
	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {

			if (i != j) {
				cout << theSynapses[j * NoTotal + i]->gsyn() << " ";
			} else
				cout << "0.0";
		}
		cout << endl;
	}

}

void neuroSynNetworkDet::enable() {
	model = new NeuronModel(&neurs, &syns, N, cerr);
	x = new double[N];
	xn = new double[N];
	machine = new rk65n(N, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);
}

void neuroSynNetworkDet::SetWeights(Array2D<double> weights,
		Array1D<double> SeqArray) {

	mWeights = weights;
	mSeqArray = SeqArray;

	for (int i = 0; i < NoTotal; i++) {
		for (int j = 0; j < NoTotal; j++) {

			if (i != j) {
				theSynapses[j * NoTotal + i]->set_gsyn(mWeights[i][j]);
			}
		}
	}
}

void neuroSynNetworkDet::ScaleWeights(double factor) {

	for (int j = 0; j < NoTotal; j++) {
		for (int i = 0; i < NoTotal; i++) {
			if (i != j)
				theSynapses[j * NoTotal + i]->set_gsyn(mWeights[i][j] * factor);

		}

	}
}

void neuroSynNetworkDet::ScaleBias(double factor) {
	double Parameters[NEUROSYN_PNO];

	for (int i = 0; i < NEUROSYN_PNO; i++)
		Parameters[i] = PARAMS_LN[i];

	for (int i = 0; i < NoTotal; i++) {
		theLNs[i]->p[11] = (theLNs[i]->p[11] * factor);
	}
}

void neuroSynNetworkDet::SetCC(Array1D<double> biasCC) {
	double Parameters[NEUROSYN_PNO];

	for (int i = 0; i < NEUROSYN_PNO; i++)
		Parameters[i] = PARAMS_LN[i];

	for (int i = 0; i < NoTotal; i++) {
		Parameters[11] = biasCC[i];
		theLNs[i]->set_p(Parameters);
	}
}

void neuroSynNetworkDet::init() {
	dt = 0.0001;
	dtx = 0.0;
	N = 0.0;

	double Initvalue[NEUROSYN_IVARNO];

	int counter = 0;
	int counter2 = 0;
	for (niter = neurs.begin(); niter != neurs.end(); niter++) {
		if (counter < NoTotal) {
			for (int i = 0; i < NEUROSYN_IVARNO; i++)
				Initvalue[i] = NEUROSYN_INIVARS[i] * RG.n() - 0.5 * RG.n();
			if (counter < NoDirectLNInput)
				Initvalue[5] = 0;
			else
				Initvalue[5] = RG.n();
			(*niter)->init(x, Initvalue);
			counter2++;

		}
		counter = counter + 1;
	}

	for (int i = 0; i < NoTotal; i++)
		directLNInput[i]->set_I(0.0);
}

Array1D<double> neuroSynNetworkDet::run(double tme, double inputCurr, Array1D<
		int> inputStart) {

	Array1D<double> endPoints(2 * NoTotal, 0.0);
	vector<double> spike_history;
	stringstream name;
	char thename[80];
	ofstream NSDataN, NSDataS, NSDataIysn, NSDataM;

	NSDataN.precision(10);
	NSDataS.precision(10);
	NSDataIysn.precision(10);
	NSDataM.precision(10);
	name.clear();


	double *tmp;

	x[0] = 0;

	double factor;




	if(doneFileCreate)
				{
				if (Manual)
					name << globalName << "CondF.dat" << ends;
				else
					name << "CondF.dat" << ends;

				name >> thename;
				NSDataN.open(thename);

				name.clear();

				if (Manual)
					name << globalName << "CondS.dat" << ends;
				else
					name << "CondS.dat" << ends;

				name >> thename;
				NSDataS.open(thename);

				name.clear();
				}

	double accumulater = 0.0;
	while (x[0] < tme) {


		if (Gradual) {

			if (int(x[0] * 10) % int(scaleBin * 10) == 0 && x[0] < scaleTime) {

				factor = (scaleBin + x[0]) / scaleTime;
				ScaleWeights(0.5 + 0.3 * factor);
		//		cout << (0.5 + 0.4 * factor) * percentCritical << endl;

			}

			if (int((x[0] - scaleStart2) * 10) % int(scaleBin * 10) == 0
					&& x[0] > scaleStart2 && x[0] < (scaleStart2 + scaleTime2)) {

				factor = (scaleBin + x[0] - scaleStart2) / scaleTime2;
				ScaleWeights(0.8  + 0.2*factor);
		//		cout << (0.9  + 0.1*factor) * percentCritical << endl;

			}
		}

		for (int i = 0; i < NoDirectLNInput; i++) {

			if (x[0] < IMPULSESTART || x[0] > IMPULSESTART + IMPULSEDUR)
				directLNInput[i]->set_I(RGaus.n() * NoiseMag);
			else
				directLNInput[i]->set_I(RGaus.n() * NoiseMag + inputCurr);
		}

		double tdt = 0.0;

		for (int i = 0; i < NoTotal; i++) {
			if (isnan(theLNs[i]->E(x))) {
				exit(1);
			}


			if (theLNs[i]->S(x) < 0.001 && x[0] > (scaleStart2 + scaleTime2)  &&  CreateFile)
				exit(1);

		//	if(accumulater > 200.0)
			//	exit(1);
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
			theLNs[i]->spike_detect(x);

		//out to screen progress
		if (int(x[0] * 10) % 5000 == 0)
			cout << x[0] << endl;

		if(doneFileCreate){

			if(x[0] > OutTime)
				NSDataN << x[0];

			for (int i = 0; i < NoTotal; i++) {
				double spiker = 0;
				if (theLNs[i]->start_spiking)
					spiker = 1;


				if(x[0] > OutTime)
					NSDataN << " " << spiker;
			}

			if(x[0] > OutTime)
				NSDataN << endl;



			if(x[0] > OutTime)
			{
				NSDataM << x[0];
			for (int i = 0; i < NoTotal; i++)
				if(x[0] > OutTime)
					NSDataM << " " << theLNs[i]->E(x);


			NSDataM << endl;
			}





			if(x[0] > OutTime)
						{

				NSDataS << x[0];
			for (int i = 0; i < NoTotal; i++) {
				NSDataS << " " << theLNs[i]->S(x);
			}
			NSDataS << endl;
						}


		}

		if (x[0] > (IMPULSESTART - 800) && x[0] < IMPULSESTART) {

			for (int i = 0; i < NoTotal; i++) {
				double spiker = 0;

				if (theLNs[i]->start_spiking)
					spiker = 1.0;

				endPoints[i] = endPoints[i] + spiker;
			}

		}

		if (x[0] > (THETIME - 800)) {

			for (int i = NoTotal; i < 2 * NoTotal; i++) {
				double spiker = 0;

				if (theLNs[i - NoTotal]->start_spiking)
					spiker = 1.0;

				endPoints[i] = endPoints[i] + spiker;
			}

		}

	}

	if (doneFileCreate) {
		NSDataN.close();
		NSDataS.close();
		NSDataM.close();
	}

	for (int i = 1; i < 2 * NoTotal; i++)
		endPoints[i] = endPoints[i] / 800.0;

	return endPoints;

}

#endif
