#pragma once

#include <iostream>
#include <vector>
#include "analysis.h"

#include "HepMC/GenEvent.h"

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
void get_final_state_particles(HepMC::GenEvent & hepMCevent, finalState& particles);


// Substructure variables
std::vector< double > SplittingScales( std::vector<fastjet::PseudoJet> const& jetVec );
std::vector< double > NSubjettiness( std::vector<fastjet::PseudoJet> const& jetVec, double const& jet_rad );
std::vector< double > NSubjettiness( std::vector<fastjet::PseudoJet> const& jetVec, double const& jet_Rmax, double const& jet_Rmin, double const& jet_Rho );

// Ghost association for large-R double b-tagging
void get_assoc_trkjets( fastjet::PseudoJet calojet, std::vector<fastjet::PseudoJet> trkjets, std::vector<fastjet::PseudoJet> &matched_trkjets, bool debug);

// Btagging
double btag_eff( double jet_pt );
double mistag_eff( double jet_pt );
double charm_eff( double jet_pt );