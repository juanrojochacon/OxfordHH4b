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
double pythiaSeed();

// Pileup
bool pileupSimulated();
int npileupEvents();
std::string minBiasFile();

// Analysis initialisation
void InitAnalyses(	std::vector<Analysis*>& HH4bAnalyses,
					std::vector<Analysis*>& signalAnalyses,
					std::vector<Analysis*>& backgroundAnalyses);
void InitSampleAnalyses( std::vector<Analysis*>& sampleAnalyses, std::string const& samplename );