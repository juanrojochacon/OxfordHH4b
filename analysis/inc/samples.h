// samples.h
#pragma once

#include <string>

#include "utils.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "Pythia8/Pythia.h"

// Initialise Pythia for reading from a LHE file, and showering it
void InitPythia(Pythia8::Pythia & pythiaRun, std::string const& eventfile, int const& nevt_max, double& weight_norm);
// Initiate a HepMC sample, and determine the cross-section normalisation
void InitHepMC( std::string const& eventfile, int const& nevt_max, double& weight_norm);

// Fetch the final state for a pythia input
void get_final_state_particles(Pythia8::Pythia & pythiaRun, finalState& particles, double& unit_weight);
// Fetch the final state particles for a HepMC file
void get_final_state_particles(std::ifstream& hepmc_is, finalState& particles, double& unit_weight);
