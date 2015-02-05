// trainingdata.h
#pragma once

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

class trainingDatum
{
public:
	trainingDatum(int const& _signal, string const& _source, double const& _weight, int const& _nKin, const double* _kinematics);
	~trainingDatum();

	// Get methods
	int const& getSignal() const {return signal;};
	string const& getSource() const {return source;};
	const double& getWeight() const {return weight;};
	const int& getNkin() const {return nKin;};
	const double* getKinematics() const {return kinematics;};

	void printDatum() const;

	void standardiseKin(const double* norm, const double* shift);

private:
	const int signal; 	 //!< Is the datum signal or background
	const string source; //!< Source sample of the datum

	const double weight; //!< Event weight

	const int nKin;		//!< Number of kinematic points
	double* kinematics; //!< Input for MVA
};

// Normalise a vector of trainingData
void norm_trainingData (vector<trainingDatum*> &dataPoints);
