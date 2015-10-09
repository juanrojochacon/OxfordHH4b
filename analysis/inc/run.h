#pragma once

#include <string>
#include <iostream>
#include <vector>

// Runcard header
class Analysis;

struct eventSample
{
	const std::string eventfile;	// Filename for event sample
	const std::string samplename;	// Identifying name for sample
	const double xsec_norm;			// xSection normalistion (e.g K-factors)
	const bool signal;				// Is the sample a signal process
	const bool hepmc;				// HepMC? (if false, use Pythia)
	const int nevt_sample;			// Number of events in sample
};

// Sample Fetching
int GetNSamples();
eventSample GetSample( int const& isample );

// Global settings
double GetPSmear();	// Return the % of Jet momenta to smear by
double GetESmear();	// Return the % of Jet energy to smear by

bool pythiaShowered(); // Is pythia showering

int& subSample(); // Current working subsample
int sampleSize(); // Size point of sample;
int sampleStart(); // Start point of sample;

// RNG seeds
double& pythiaSeed();
double& systemSeed();

// Pileup
bool pileupSimulated();
int npileupEvents();
int npileupTotal();
std::string minBiasFile();

void InitSampleAnalyses( std::vector<Analysis*>& sampleAnalyses, std::string const& samplename, int const& subsample );