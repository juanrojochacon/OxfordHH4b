#pragma once

#include <iostream>
#include <vector>
#include "analysis.h"

namespace Pythia8{
class Pythia;
}

/*
This routine initializes Pythia8
with all the settings for the shower and underlying event
 */
void InitPythia(Pythia8::Pythia & pythiaRun, string eventfile);

/*
This routine initializes Pythia8 with no shower/UE/Hadronisation
 */
void InitPythia_PartonLevel(Pythia8::Pythia & pythiaRun, string eventfile);


/*
  Get the information on all final state particles
 */
void get_final_state_particles(Pythia8::Pythia & pythiaRun, finalState& particles);
