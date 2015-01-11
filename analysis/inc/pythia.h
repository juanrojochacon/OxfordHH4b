#pragma once

#include <string>

#include "utils.h"
#include "settings.h"

namespace Pythia8{
class Pythia;
}

/*
This routine initializes Pythia8
with all the settings for the shower and underlying event
 */
void InitPythia(Pythia8::Pythia & pythiaRun, std::string const& eventfile, int const& nevt_max, double& weight_norm);

// Fetch the final state for a pythia input
void get_final_state_particles(Pythia8::Pythia & pythiaRun, finalState& particles, double& unit_weight);
