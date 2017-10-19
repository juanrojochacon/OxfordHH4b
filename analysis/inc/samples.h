// samples.h
#pragma once

#include <string>
#include <vector>

#include "run.h"
#include "utils.h"

#include "HepMC/GenCrossSection.h"
#include "HepMC/GenEvent.h"
#include "Pythia8/Pythia.h"

// Initialise Pythia for reading from a LHE file, and showering it
void InitPythia(runCard const& rc, sampleCard const& sc, uint32_t const& seed,
                Pythia8::Pythia& pythiaRun, double& weight_norm);
// Initiate a HepMC sample, and determine the cross-section normalisation
void InitHepMC(runCard const& rc, sampleCard const& sc, double& weight_norm,
	       int& index_subsample_size, std::vector<long long>& subsample_indices);

// Fetch the final state for a pythia input
void get_final_state_particles(Pythia8::Pythia& pythiaRun, finalState& particles,
                               double& unit_weight);
// Fetch the final state particles for a HepMC file
void get_final_state_particles(std::ifstream& hepmc_is, finalState& particles, double& unit_weight, bool print = false);
