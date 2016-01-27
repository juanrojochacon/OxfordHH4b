// samples.h
#pragma once

#include <string>

#include "utils.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "Pythia8/Pythia.h"

// Class to parse event sample cards
class sampleCard
{
public:
	sampleCard(std::string const& filename); // Constructor from filename

	const std::string eventfile;	// Filename for event sample
	const std::string samplename;	// Identifying name for sample
	const std::string format;		// HepMC/LHE
	const double xsec_norm;			// xSection normalistion (e.g K-factors)
	const double sqrts;				// Centre of mass energy
	const bool signal;				// Is the sample a signal process
	const int nevt_sample;			// Number of events in sample
};

// Initialise Pythia for reading from a LHE file, and showering it
void InitPythia(Pythia8::Pythia & pythiaRun, std::string const& eventfile, int const& nevt_max, double& weight_norm);
// Initiate a HepMC sample, and determine the cross-section normalisation
void InitHepMC( std::string const& eventfile, int const& nevt_max, double& weight_norm);

// Fetch the final state for a pythia input
void get_final_state_particles(Pythia8::Pythia & pythiaRun, finalState& particles, double& unit_weight);
// Fetch the final state particles for a HepMC file
void get_final_state_particles(std::ifstream& hepmc_is, finalState& particles, double& unit_weight);
