/////////////////////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

#include <algorithm>
#include <cmath>

#include "parametrisation.h"
#include "trainingdata.h"
#include "random.h"


using std::cout;
using std::endl;



////////////////////////////////////////////////////////////////////////////////////////

int main() 
{  
	// Initialise RNG
	rng_init(235243356345);

	// Read datafile
	string dataPath = "./"+ string(DATADIR) + "/trainingData.dat";
	cout << "Reading data from "<< dataPath<<endl;

	ifstream datafile(dataPath.c_str());
	string line; stringstream linestream;
	getline(datafile,line);
	linestream.str(line);

	// Dummy
	for (int i=0; i<3; i++)
		linestream >> line;

	// Count entries
	int nKin = 0;
	while (linestream >> line)
		nKin++;
	
	cout <<nKin << " kinematic points found"<< endl;

	// Read entries
	vector<trainingDatum*> trainingData;
	int sigCount = 0;
	int bkgCount = 0;
	while (getline(datafile, line))
	{
		bool signal;
		string source;
		double* kinematics = new double[nKin];

		stringstream datstream(line);
		datstream >> signal >> source;
		for (int i=0; i<nKin; i++)
			datstream >> kinematics[i];

		// Add datapoint to vector
		trainingData.push_back(new trainingDatum(signal,source,nKin,kinematics));

		if (signal)
			sigCount++;
		else
			bkgCount++;

		delete[] kinematics;
	}

	cout << trainingData.size() << " datapoints found in training set."<<endl; 
	cout << sigCount<<" signal points, " <<bkgCount << " background points." <<endl;

	const double sig_wgt = ( (double) bkgCount )/( (double) sigCount );
	const double bkg_wgt = ( (double) sigCount )/( (double) bkgCount );
	cout << "Background weight: "<< bkg_wgt << " signal weight: " <<sig_wgt <<endl;

	norm_trainingData(trainingData);
	norm_trainingData(trainingData);

	// Initialise MLP
	std::vector<int>  nnArch;
	nnArch.push_back(nKin);
	nnArch.push_back(5);
	nnArch.push_back(3);
	nnArch.push_back(1);

	// Fit NN
	MultiLayerPerceptron mlp(nnArch);

	// Compute initial probabilities
	cout << "******************************************************"<<endl;
	double* outProb = new double[1];
	double fitness = 0;
	for (size_t i=0; i<trainingData.size(); i++)
	{
		*outProb = 0;
		mlp.Compute(trainingData[i]->getKinematics(), outProb);

		// Compute cross-entropy
		const double t = trainingData[i]->getSignal();
		const double tpr = *outProb;
		const double wgt = 1.0;//trainingData[i]->getSignal() ? sig_wgt:bkg_wgt;

		//fitness -= wgt*t*log(tpr)+(1.0-t)*log(1.0-tpr); // cross-entropy
		fitness += wgt*pow((t-tpr),2.0);	// MSE

	}

	cout << "Init EC: " << fitness<<endl;
	cout << "******************************************************"<<endl;

	//const int nGen = 500000;
	const int nGen = 50000;
	for (int i=0; i< nGen; i++)
	{
		MultiLayerPerceptron mutant(mlp);
		// Mutate
		const int NParam =  mutant.GetNParameters();
        mutant.GetParameters()[rng_uniform(NParam)]+=rng_gaussian(0.2);

        // Compute EC
        double mut_fitness = 0;
		for (size_t i=0; i<trainingData.size(); i++)
		{
			*outProb = 0;
			mutant.Compute(trainingData[i]->getKinematics(), outProb);

			// Compute cross-entropy
			const double t = trainingData[i]->getSignal();
			const double tpr = *outProb;
			const double wgt = 1.0;//trainingData[i]->getSignal() ? sig_wgt:bkg_wgt;

			//mut_fitness -= wgt*t*log(tpr)+(1.0-t)*log(1.0-tpr); // cross-entropy
			mut_fitness += wgt*pow((t-tpr),2.0);	// MSE
		}

		if (mut_fitness < fitness)
		{
			fitness = mut_fitness;
			mlp.CopyPars(&mutant);
		}

		if (i% 5000 == 0)
		cout << i<<"\t"<<fitness<<endl;
	}

	// Export for analysis
	stringstream arch;
	for (int i=0; i<nnArch.size(); i++)
	{
		arch << nnArch[i];
		if (i != nnArch.size()-1)
			arch << "X";
	}
	arch << "_"<<nGen<<"-Gen";

	const string mvafile = "./" + string(RESDIR) + "/nn_" + arch.str() + "_MSE.dat";
	ofstream mvaout(mvafile.c_str());

	cout << "******************************************************"<<endl;
	for (size_t i=0; i<trainingData.size(); i++)
	{
		*outProb = 0;
		mlp.Compute(trainingData[i]->getKinematics(), outProb);
		mvaout << trainingData[i]->getSource()<<"\t"<<trainingData[i]->getSignal() <<"\t"<<*outProb<<endl;
	}

	mvaout.close();

	cout << "******************************************************"<<endl;

	// End of the main progream
	return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
