#pragma once

#include <iostream>
#include <vector>
#include "analysis.h"

#include "HepMC/GenEvent.h"
#include "fastjet/PseudoJet.hh"

// Helpful typedef
typedef enum {NTAG, BTAG, CTAG, LTAG} btagType;

// Analysis utility functions
double getDPhi(double phi1, double phi2);

// Ghost association for large-R double b-tagging
void get_assoc_trkjets( fastjet::PseudoJet calojet, std::vector<fastjet::PseudoJet> trkjets, std::vector<fastjet::PseudoJet> &matched_trkjets, bool debug);

// Btagging
double btag_eff( double jet_pt );
double mistag_eff( double jet_pt );
double charm_eff( double jet_pt );

// Substructure variables
std::vector< double > SplittingScales( std::vector<fastjet::PseudoJet> const& jetVec );
double SplittingScales( fastjet::PseudoJet const& jet );

// Angular variables
double Chi(fastjet::PseudoJet h0, fastjet::PseudoJet h1);
