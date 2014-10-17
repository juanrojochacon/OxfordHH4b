/*
This are the common settings used in the analysis code
 */

// Cuts for the b-jet candidates for b-tagging
double const pt_btagging=15;

// hadronic dijet mass resolution -> Taken to be 15 GeV
// Taken from LHC H -> bb papers
double const mass_resolution = 0.12; // 15 GeV 

// Higgs mass
double const m_higgs = 125.0;

// Jet definition
// Small radius jet reconstruction
double const jetR_0p5=0.5; // To avoid overlapping b's as much as possible
// Large radius jet reconstruction
double const jetR_1p2=1.2; // To try to merge two b quarks into the same jet

JetDefinition CA10(cambridge_algorithm, jetR_1p2);

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

// To check energy-momentum conservation
double const Eref=14000; // Samples generated for LHC 14 TeV
double const tol_emom=1.0;

//-------------------------------------------------------------

