#pragma once

#include <iostream>
#include <vector>
#include "analysis.h"

#include "HepMC/GenEvent.h"

// Analysis utility functions

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