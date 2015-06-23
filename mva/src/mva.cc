/////////////////////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

#include <algorithm>
#include <cmath>
#include <limits>


#include "parametrisation.h"
#include "trainingdata.h"
#include "random.h"

#include <signal.h>

using std::cout;
using std::endl;

bool interrupt = false;

void catch_int( int signum )
{
	interrupt = true;
}

////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{  

  	signal(SIGINT, catch_int);

	// Initialise RNG
	rng_init(235243356345);

	if (argc != 2)
	{
		cerr << "Error: Wrong number of arguments!"<<endl;
		cerr << "Usage: mva <path_to_trainingdata>" <<endl;
		exit(-1);
	}

	// Read datafile
	const string dataPath = argv[1];
	const string dataName = dataPath.substr(2,dataPath.length()-6);
	cout << "Reading data from "<< dataPath<< "  Named: " << dataName<< endl;

	ifstream datafile(dataPath.c_str());
	string line; stringstream linestream;
	getline(datafile,line);
	linestream.str(line);

	// Dummy indices
	for (int i=0; i<4; i++)
		linestream >> line;

	// Count entries
	int nKin = 0;
	stringstream kinstream;
	while (linestream >> line)
	{
		kinstream << line<<"\t";
		nKin++;
	}

	cout <<nKin << " kinematic points found"<< endl;

	// Read entries
	vector<trainingDatum*> trainingData;
	int sigCount = 0;
	int bkgCount = 0;

	double sigWeight = 0;
	double bkgWeight = 0;

	while (getline(datafile, line))
	{
		bool signal;
		string source;
		double weight;

		double* kinematics = new double[nKin];

		stringstream datstream(line);
		datstream >> signal >> source >> weight;
		for (int i=0; i<nKin; i++)
			datstream >> kinematics[i];

		// Add datapoint to vector
		trainingData.push_back(new trainingDatum(signal,source,weight,nKin,kinematics));

		if (signal)
		{
			sigCount++;
			sigWeight+=weight;
		}
		else
		{
			bkgCount++;
			bkgWeight+=weight;
		}

		delete[] kinematics;
	}

	cout << trainingData.size() << " datapoints found in training set."<<endl; 
	cout << sigCount<<" signal points, " <<bkgCount << " background points." <<endl;
	cout << sigWeight<<" signal weight, " <<bkgWeight << " background weight." <<endl;

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
	double fitness = std::numeric_limits<double>::infinity();

	const int nGen = 50000;
	for (int i=0; i< nGen; i++)
	{
		if (interrupt) break;

		// Mutate
		MultiLayerPerceptron mutant(mlp);
		const int NParam =  mutant.GetNParameters();
		mutant.GetParameters()[rng_uniform(NParam)]+=rng_gaussian(0.2);

		// Compute mutant fitness
		double mut_fitness = 0;
		for (size_t j=0; j<trainingData.size(); j++)
		{
			*outProb = 0;
			mutant.Compute(trainingData[j]->getKinematics(), outProb);

			// Compute cross-entropy
			const double t = trainingData[j]->getSignal();
			const double tpr = *outProb;
			const double wgt = trainingData[j]->getWeight() / (trainingData[j]->getSignal() ? (sigWeight/sigCount):(bkgWeight/bkgCount));

			mut_fitness -= wgt*t*log(tpr)+(1.0-t)*log(1.0-tpr); // cross-entropy
			//mut_fitness += wgt*pow((t-tpr),2.0);	// MSE

			if (mut_fitness > fitness) break;
		}

		// Selection
		if (mut_fitness < fitness)
		{
			fitness = mut_fitness;
			mlp.CopyPars(&mutant);
		}

		// Write progress to screen
		if (i% 50 == 0)
			cout << i<<"\t"<<fitness<<endl;
	}

	// Export for analysis
	stringstream arch;
	for (size_t i=0; i<nnArch.size(); i++)
	{
		arch << nnArch[i];
		if (i != nnArch.size()-1)
		arch << "X";
	}
	arch << "_"<<nGen<<"-Gen";

	const string mvafile = "./" + string(RESDIR) + "/nn_" + arch.str() + "_"+ dataName+ ".dat";
	ofstream mvaout(mvafile.c_str());

	cout << "******************************************************"<<endl;
	for (size_t i=0; i<trainingData.size(); i++)
	{
		*outProb = 0;
		mlp.Compute(trainingData[i]->getKinematics(), outProb);
		mvaout << trainingData[i]->getSource()<<"\t"<<trainingData[i]->getSignal() <<"\t"<<trainingData[i]->getWeight()<<"\t"<<*outProb<<endl;
	}

	mvaout.close();
	cout << "******************************************************"<<endl;

	const string netfile = "./" + string(RESDIR) + "/nn_" + arch.str()+ "_"+ dataName+  ".net";
	mlp.Export(netfile, kinstream.str());

	// End of the main progream
	return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
