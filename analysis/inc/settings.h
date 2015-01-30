/*
Common settings used in the analysis code
 */

#pragma once


// FastJet
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
// FastJet contrib
//#include "fastjet/contrib/VariableRPlugin.hh"

// hadronic dijet mass resolution -> Taken to be 15 GeV
// Taken from LHC H -> bb papers
double const mass_resolution = 0.12; // 15 GeV 

// Higgs mass
double const m_higgs = 125.0;

// b tagging
// Choose working point with high purity
double const btag_prob = 0.80; // Probability of correct b tagging
double const btag_mistag = 0.01; // Mistag probability  
const double ctag_prob = 0.17; // c-mistag rate ~1/6.

// Run II Lumi in 1/fb
double const lumi_run2=300;
double const lumi_hllhc=3000;

//-------------------------------------------------------------
