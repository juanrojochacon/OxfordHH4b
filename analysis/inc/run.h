#pragma once

#include <iostream>
#include <string>
#include <vector>

// Runcard header - parsing of runcards and global settings
class Analysis;

// Class to parse event sample cards
class sampleCard {
  public:
    sampleCard(std::string const& filename); // Constructor from filename
    const std::string eventfile;             // Filename for event sample
    const std::string eventpath;             // Path to event sample
    const std::string samplename;            // Identifying name for sample
    const std::string format;                // HepMC/LHE
    const bool        hepmc;                 // isHepMC
    const double      xsec_norm;             // xSection normalistion (e.g K-factors)
    const double      sqrts;                 // Centre of mass energy
    const bool        is_signal;             // Is the sample a signal process
    const int         nevt_sample;           // Number of events in sample
};

class runCard {
  public:
    runCard(std::string const& filename);
    const std::string runname;        // Run name
    const int         sub_samplesize; // Size of subsample
    const int         npileup;        // Number of pileup events
    const double      jetEsmear;      // Jet energy smearing (%)
    const bool        pythiaShower;   // Shower LHE events
    const int         runseed;        // Global seed
};

// Initialise list of analyses
void InitAnalyses(std::vector<Analysis*>& sampleAnalyses, runCard const& run,
                  sampleCard const& sample, int const& subsample);
