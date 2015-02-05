// trainingdata.cc

#include "trainingdata.h"
#include <algorithm>

trainingDatum::trainingDatum(int const& _signal, string const& _source, double const& _weight, int const& _nKin, const double* _kinematics):
signal(_signal),
source(_source),
weight(_weight),
nKin(_nKin)
{
	kinematics = new double[nKin];
	for (int i=0; i<nKin; i++)
		kinematics[i] = _kinematics[i];
}

trainingDatum::~trainingDatum()
{
	delete[] kinematics;
}


void trainingDatum::printDatum() const
{
	cout << signal << "\t" << source;
	for (int i=0; i<nKin; i++)
		cout << "\t" <<kinematics[i];
	cout <<endl;
}

void trainingDatum::standardiseKin(const double* norm, const double* shift)
{
	for (int i=0; i<nKin; i++)
		kinematics[i] = norm[i]*kinematics[i] + shift[i];
}


void norm_trainingData (vector<trainingDatum*> &dataPoints)
{
	if (dataPoints.size() == 0)
	{
		cerr << "norm_trainingData Error: no datapoints!"<<endl;
		exit(-1);
	}

	const int nKin = dataPoints[0]->getNkin();

	double* kinNorm  = new double[nKin];
	double* kinShift = new double[nKin];

	for (int i=0; i<nKin; i++)
	{
		vector<double> kinVals;
		for (size_t j=0; j<dataPoints.size(); j++)
			kinVals.push_back(dataPoints[j]->getKinematics()[i]);

		const double maxval = *std::max_element(kinVals.begin(), kinVals.end());
		const double minval = *std::min_element(kinVals.begin(), kinVals.end());

		// Normalise to \pm sqrt(3)
		kinNorm[i] = 2.0*sqrt(3)/(maxval - minval);
		kinShift[i] =  - ( kinNorm[i]*minval + sqrt(3) );
	}

	for (size_t i=0; i<dataPoints.size(); i++)
		dataPoints[i]->standardiseKin(kinNorm, kinShift);

}
