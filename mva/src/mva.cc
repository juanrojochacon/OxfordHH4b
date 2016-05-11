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
const double trsplit = 0.5;

void catch_int( int signum )
{
	interrupt = true;
}

void ComputeFitness(MultiLayerPerceptron const& mlp, trainingDatum const* datum, double const& weight, double& fitness)
{
}

////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{  

  	signal(SIGINT, catch_int);

	if (argc != 3)
	{
		cerr << "Error: Wrong number of arguments!"<<endl;
		cerr << "Usage: mva <path_to_trainingdata> <boostrap>" <<endl;
		exit(-1);
	}

	// Read datafile
	const string dataPath = argv[1];
	const string dataName = dataPath.substr(2,dataPath.length()-6);
	const string bootStrap = argv[2];
	cout << "Reading data from "<< dataPath<< "  Named: " << dataName<< endl;

	// Initialise RNG
	rng_init(time(NULL)/atoi(bootStrap.c_str()));

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

	int sigCount = 0;
	int bkgCount = 0;

	double sigWeight = 0;
	double bkgWeight = 0;

	vector<trainingDatum*> totalData;
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

		totalData.push_back(new trainingDatum(signal,source,weight,nKin,kinematics));

		delete[] kinematics;
	}

	cout << totalData.size() << " datapoints found in total."<<endl; 
	cout << sigCount<<" signal points, " <<bkgCount << " background points." <<endl;
	cout << sigWeight<<" signal weight, " <<bkgWeight << " background weight." <<endl;

	// Normalise
	norm_trainingData(totalData);

	vector<trainingDatum*> trainingData;
	vector<trainingDatum*> validationData;

	for (size_t i=0; i<totalData.size(); i++)
	{
		if ( rng_uniform() < trsplit)
			trainingData.push_back(totalData[i]);
		else
			validationData.push_back(totalData[i]);
	}

	cout << "Cross-validation: "<<endl;
	cout << trainingData.size() << " datapoints found in training set."<<endl; 
	cout << validationData.size() << " datapoints found in validation set."<<endl; 


	// Initialise MLP
	std::vector<int>  nnArch;
	nnArch.push_back(nKin);
	nnArch.push_back(5);
	nnArch.push_back(3);
	nnArch.push_back(1);

	MultiLayerPerceptron mlp(nnArch); 	 // Fit NN
	MultiLayerPerceptron cv_mlp(nnArch); // Cross-validated NN

	// Compute initial probabilities
	cout << "******************************************************"<<endl;
	double* outProb = new double[1];
	double fitness = std::numeric_limits<double>::infinity();
	double lookback_fitness = std::numeric_limits<double>::infinity();
	int lookback_gen = 0;

	const int nGen = 50000;
	for (int i=0; i< nGen; i++)
	{
		if (interrupt) break;

		// Mutate
		MultiLayerPerceptron mutant(mlp);
		const int NParam =  mutant.GetNParameters();
		mutant.GetParameters()[rng_uniform(NParam)]+=rng_gaussian(10.0)/pow(i+1, 0.5);

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
			if (mut_fitness > fitness) break;
		}

		// Compute mutant CV
		double mut_lookback_fitness = 0;
		if (mut_fitness < fitness)
			for (size_t j=0; j<validationData.size(); j++)
			{
				*outProb = 0;
				mutant.Compute(validationData[j]->getKinematics(), outProb);

				// Compute cross-entropy
				const double t = validationData[j]->getSignal();
				const double tpr = *outProb;
				const double wgt = validationData[j]->getWeight() / (validationData[j]->getSignal() ? (sigWeight/sigCount):(bkgWeight/bkgCount));

				mut_lookback_fitness -= wgt*t*log(tpr)+(1.0-t)*log(1.0-tpr); // cross-entropy
				if (mut_lookback_fitness > lookback_fitness) break;
			}

		// Selection
		if (mut_fitness < fitness)
		{
			fitness = mut_fitness;
			mlp.CopyPars(&mutant);

			if (mut_lookback_fitness < lookback_fitness)
			{
				lookback_fitness = mut_lookback_fitness;
				lookback_gen = i;
				cv_mlp.CopyPars(&mutant);
			}
		}

		// Write progress to screen
		if (mut_fitness == fitness && (i % 10) == 0 )
		{
			const double reldiff = abs(mut_lookback_fitness-lookback_fitness)/lookback_fitness;
			cout << i<<"\t"<<fitness<<"\t"<<mut_lookback_fitness<<"  "<<lookback_fitness<<"  "<<reldiff<<endl;
		}
	}

	cout << "Completed, "<<nGen<<" generations: stopped at :"<<lookback_gen<<endl;

	// Export for analysis
	stringstream arch;
	for (size_t i=0; i<nnArch.size(); i++)
	{
		arch << nnArch[i];
		if (i != nnArch.size()-1)
		arch << "X";
	}
	arch << "_"<<nGen<<"-Gen";

	const string mvafile = "./" + string(RESDIR) + "/nn_" + arch.str() + "_"+ dataName+"_"+bootStrap+ ".dat";
	cout << "Exporting to: "<<mvafile<<endl;
	ofstream mvaout(mvafile.c_str());

	cout << "******************************************************"<<endl;
	for (size_t i=0; i<totalData.size(); i++)
	{
		*outProb = 0;
		cv_mlp.Compute(totalData[i]->getKinematics(), outProb);
		mvaout << totalData[i]->getSource()<<"\t"<<totalData[i]->getSignal() <<"\t"<<totalData[i]->getWeight()<<"\t"<<*outProb<<endl;
	}

	mvaout.close();
	cout << "******************************************************"<<endl;

	const string netfile = "./" + string(RESDIR) + "/nn_" + arch.str()+ "_"+ dataName+"_"+bootStrap+  ".net";
	const string parfile = "./" + string(RESDIR) + "/nn_" + arch.str()+ "_"+ dataName+"_"+bootStrap+   ".par";
	cv_mlp.Export(netfile, kinstream.str());
	cv_mlp.ExportPars(parfile);

	// End of the main progream
	return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
