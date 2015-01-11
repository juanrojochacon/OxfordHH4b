#pragma once

#include <string>

#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"

#include "utils.h"
#include "settings.h"


// Initiate a HepMC sample, and determine the cross-section normalisation
void InitHepMC( std::string const& eventfile, int const& nevt_max, double& weight_norm);

// Fetch the final state particles for a HepMC file
void get_final_state_particles(std::ifstream& hepmc_is, finalState& particles, double& unit_weight);
