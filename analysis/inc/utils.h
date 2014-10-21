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
  Get the information on all final state particles
 */
void get_final_state_particles(Pythia8::Pythia & pythiaRun, finalState& particles);

/*
This routine simulates the b tagging
In its present form, it only requires that at least
one of the jet constitutents is a b jet
b quarks must be above some minimum pt
 */
bool btagging(int ijet, std::vector<fastjet::PseudoJet> jets);