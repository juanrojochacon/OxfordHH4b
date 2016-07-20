/////////////////////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

#include <algorithm>
#include <cmath>
#include <limits>

#include "parametrisation.h"
#include "trainingdata.h"

#include <signal.h>

using std::cout;
using std::endl;

////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{  
	if (argc != 4)
	{
		cerr << "Error: Wrong number of arguments!"<<endl;
		cerr << "Usage: apply_mva <path_to_trainingdata> <path to mva> <ANNcut>" <<endl;
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

		totalData.push_back(new trainingDatum(signal,source,weight,nKin,kinematics));

		delete[] kinematics;
	}

	// Normalise
	norm_trainingData(totalData);
	cout << totalData.size() << " datapoints found in total."<<endl; 

	// Initialise MLP
	std::vector<int>  nnArch;
	nnArch.push_back(nKin);
	nnArch.push_back(5);
	nnArch.push_back(3);
	nnArch.push_back(1);

	// Read NN
	const string parPath = argv[2];
	cout << "Reading data from "<< parPath<< endl;
	MultiLayerPerceptron mlp(nnArch); 	 // Fit NN
	mlp.ImportPars(parPath);


	// Read cut point
	const double ANNcut = atof(argv[3]);
	cout << "ANN cut at: "<< ANNcut <<std::endl;

	cout << "******************************************************"<<endl;
	double *outProb = new double;
	double passWeight = 0;
	int npass = 0;
	for (size_t i=0; i<totalData.size(); i++)
	{
		*outProb = 0;
		mlp.Compute(totalData[i]->getKinematics(), outProb);
		if (*outProb > ANNcut)
		{
			passWeight += totalData[i]->getWeight();
			npass++;
		}
	}
	cout << "******************************************************"<<endl;

	std::ofstream res("results.dat");
	res << std::scientific << std::setprecision(6);
	res << "4.000000e+00    5.000000e+00    " << passWeight <<"    "<<  passWeight/sqrt(npass) <<"    "<< passWeight/sqrt(npass)<<std::endl; 
	res.close();
	// End of the main progream
	return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
