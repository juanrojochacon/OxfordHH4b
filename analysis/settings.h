/*
This are the common settings used in the analysis code
 */

// b jets
double const pt_bjet = 25;
double const eta_bjet=2.5;

// Cuts for the b-jet candidates for b-tagging
double const pt_btagging=20;

// Minimum pt for the hh system
double const pthh_cut=40;

// hadronic dijet mass resolution -> Taken to be 15 GeV
// Taken from LHC H -> bb papers
double const mass_resolution = 0.12; // 15 GeV 

// Higgs mass
double const m_higgs = 125.0;

// Jet definition
double const jetR=0.5; // To avoid overlapping b's as much as possible

// Jet definition for substructure
double const Rsb = 1.2;
JetDefinition CA10(cambridge_algorithm, Rsb);

// Set parameters of the mass drop tagger for jet substructure
// mu = 0.67 and y = 0.09 are the default choice in FastJet
double const mu = 0.67;
double const ycut = 0.09;
// The choice ycut = 0.15 is also recommended by Gavin to optimize S/sqrt(B)

// b tagging
// Choose working point with high purity
double const btag_prob = 0.80; // Probability of correct b tagging
double const btag_mistag = 0.01; // Mistag probability  

// Run II Lumi in 1/fb
double const lumi_run2=300;
double const lumi_hllhc=3000;

//-------------------------------------------------------------

