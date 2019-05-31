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

	NeuroSyn **theLNs;
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
	fullMGC();
	~fullMGC();

	void generateNetwork();
	void generateConnect();
	void init();
	void enable();
	Array1D<double> run(double tme, double inputCurr,Array1D<int> inputStart);
	double TransferInv(double in, double *params);
	void SetWeights(Array2D<double> weights, Array1D<double> SeqArray);
	void SetCC(Array1D<double> biasCC);

};

fullMGC::fullMGC() {
	dt = 0.0001;
	dtx = 0.0;
	N = 0.0;

	directLNInput = new DCInput*[NoTotal];
	theLNs = new NeuroSyn*[NoTotal];

	directORNInput = new DCInput*[NoORNs];
	theORNs = new NeuroSynRate*[NoORNs];

	//Set up network memory and space


	for (int i = 0; i < NoTotal; i++) {
		theLNs[i] = new NeuroSyn(i, PARAMS_LN);
		neurs.push_back(theLNs[i]);
		directLNInput[i] = new DCInput(theLNs[i], 1.0);
	}




	//set up connections
		theSynapses = new Emptysynapse*[NoTotal*NoTotal];
		for (int i = 0; i < NoTotal; i++) {
					for (int j = 0; j < NoTotal; j++) {
						if (i != j) {
							theSynapses[j*NoTotal+i] = new Emptysynapse(theLNs[i], theLNs[j],0.0);
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
					aSynapse = new Emptysynapse(theORNs[i],theLNs[i], -ORN2LN);

				}

	//enable integaror
	enable();
}

fullMGC::~fullMGC() {
	list<neuron *>::iterator i;
	list<synapse *>::iterator j;
	for (i = neurs.begin(); i != neurs.end(); i++) {
		delete *i;
	}
	for (j = syns.begin(); j != syns.end(); j++) {
	}

	delete[] theLNs;
		delete[] theORNs;
		delete[] directLNInput;
		delete[] directORNInput;
}

void fullMGC::enable() {
	model = new NeuronModel(&neurs, &syns, N, cerr);
	x = new double[N];
	xn = new double[N];
	machine = new rk65n(N, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);
}



void fullMGC::SetWeights(Array2D<double> weights, Array1D<double> SeqArray)
{


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


void fullMGC::SetCC(Array1D<double> biasCC)
{
	double Parameters[NEUROSYN_PNO];

		for (int i = 0; i < NEUROSYN_PNO; i++)
				Parameters[i] = PARAMS_LN[i];


		for (int i = 0; i < NoTotal; i++) {
			Parameters[11] = biasCC[i];
			theLNs[i]->set_p(Parameters);
		}
}


void fullMGC::init() {
	dt = 0.0001;
	dtx = 0.0;
	N = 0.0;

	double Initvalue[NEUROSYN_IVARNO];


	int counter = 0;
	int counter2 = 0;
	for (niter = neurs.begin(); niter != neurs.end(); niter++) {
		if (counter <  NoTotal+NoPNs) {
			for (int i = 0; i < NEUROSYN_IVARNO; i++)
					Initvalue[i] = NEUROSYN_INIVARS[i]*RG.n()-0.5*RG.n();
			if(counter< NoDirectLNInput)
			Initvalue[5] = 0;
			else
				Initvalue[5] = RG.n();				//mSeqArray[counter2];
			(*niter)->init(x, Initvalue);
			counter2++;

		}
		else{
			double temp[1];
			temp[0] = 0;
			(*niter)->init(x, temp);
		}

		counter = counter + 1;
	}

	for (int i = 0; i < NoTotal; i++)
			directLNInput[i]->set_I(0.0);
	for (int i = 0; i < NoORNs; i++)
			directORNInput[i]->set_I(0.0);
}

Array1D<double> fullMGC::run(double tme, double inputCurr,Array1D<int> inputStart) {



	Array1D<double> endPoints(2*NoTotal, 0.0);
	vector<double> spike_history;
	stringstream name;
	char thename[80];
	ofstream NSDataN, NSDataS, NSDataM;

	NSDataN.precision(10);
	NSDataS.precision(10);
	NSDataM.precision(10);
	name.clear();




	if(plotInc)
			name << globalName << "DataF" << int(10000000*inputCurr) << ".dat"<< ends;
		else
		{
			if(Manual)
				name << globalName   <<  "CondF"<< "N"  << "fs"  <<  ".dat" << ends;
			else
				 name << "CondF"<< "N"  << "fs"  <<  ".dat" << ends;

		}


	name >> thename;
		if (OutDat)
			NSDataN.open(thename);

		name.clear();

		if(plotInc)
		name << globalName   << "DataS" << int(10000000*inputCurr) << ".dat"<< ends;
		else{

			if(Manual)
				name << globalName <<  "CondS"<< "N"  << "fs"  <<  ".dat" << ends;
			else
				name << "CondS"<< "N"  << "fs"  <<  ".dat" << ends;
		}

		name >> thename;
		if (OutDat)
			NSDataS.open(thename);

		name.clear();

		if(plotInc)
		name << globalName << "DataM" << int(10000000*inputCurr) << ".dat"<< ends;
		else
		{
			if(Manual)
				name << globalName  << "CondM"<< "N"  << "fs"  <<  ".dat" << ends;
			else
				name << "CondM"<< "N"  << "fs"  <<  ".dat" << ends;

		}
	name >> thename;
	if (OutDat)
		NSDataM.open(thename);

	double *tmp;

	x[0] = 0;
	while (x[0] < tme) {

		for (int i = 0; i < NoTotal; i++) {
			//		cerr << theLNs[i]->E(x) << " ";
			if (isnan(theLNs[i]->E(x))) {
				//		cerr << "nan encountered!" << endl;
				exit(1);
			}
		}
		//		cerr << endl;


		for (int i = 0; i < 7; i++) {
			cout << x[i] << " ";
		}


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

		for (int i = 0; i < NoTotal; i++)
			theLNs[i]->spike_detect(x);



		if(!plotInc)
			cout << x[0] << endl;

		if (OutDat) {

			//	cout << x[0] << endl;
			NSDataN << x[0];



			for (int i = 0; i < NoTotal; i++) {
				double spiker = 0;
				if (theLNs[i]->start_spiking)
					spiker = 1;

				NSDataN << " " << spiker;
			}





			NSDataN << endl;

			NSDataM << x[0];


			for (int i = 0; i < NoTotal; i++)
				NSDataM << " " << theLNs[i]->E(x);





			NSDataM << endl;

			NSDataS << x[0];


			for (int i = 0; i < NoTotal; i++) {
				NSDataS << " " << theLNs[i]->S(x);
			}





			NSDataS << endl;
		}



		if (x[0] > (IMPULSESTART - 800) && x[0] < IMPULSESTART) {

			for (int i = 0; i < NoTotal; i++) {
				double spiker = 0;
		//		theLNs[i]->spike_detect(x);
				if (theLNs[i]->start_spiking)
					spiker = 1.0;

				endPoints[i] = endPoints[i] + spiker;
			}

		}


		if (x[0] > (THETIME - 800)) {

			for (int i = NoTotal; i < 2*NoTotal; i++) {
				double spiker = 0;
	//			theLNs[i-NoTotal]->spike_detect(x);
				if (theLNs[i-NoTotal]->start_spiking)
					spiker = 1.0;

				endPoints[i] = endPoints[i] + spiker;
			}

		}

	}

	if (OutDat) {
		NSDataN.close();
		NSDataS.close();
		NSDataM.close();
	}

	for (int i = 1; i < 2*NoTotal; i++)
		endPoints[i] = endPoints[i]/800.0;

	return endPoints;

}

double fullMGC::TransferInv(double in, double *params) {

	double alpha = params[14];
		double beta = params[15];
		double tr = params[16];

		double lF = beta/log((alpha*exp(beta*tr))/in-alpha/in+1);

		double out =  (lF-C)/M;
		return out;

}
#endif
