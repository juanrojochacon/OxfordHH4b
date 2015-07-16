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

// Mass-drop tagger
double const mu = 0.67;
double const ycut = 0.09;

// b tagging
// Choose working point with high purity
double const btag_prob = 0.80; // Probability of correct b tagging
double const btag_mistag = 0.01; // Mistag probability  
double const ctag_prob = 0.01;//0.17; // c-mistag rate ~1/6.

// bb tagging
// Choose working point with high purity
double const bbtag_prob = 0.80; // Probability of correct bb tagging
double const bbtag_mistag = 0.01; // Mistag probability

// Run II Lumi in 1/fb
double const lumi_run2=300;
double const lumi_hllhc=3000;

// Initial seed for smearing
uint const smear_seed = 23478234;

typedef enum {NTAG, BTAG, CTAG, LTAG} btagType;
//-------------------------------------------------------------
