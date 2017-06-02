// detector.h
#pragma once

#include "analysis.h"
#include "pcg/pcg_random.hpp"
#include "run.h"

#include "Pythia8/Pythia.h"

// Run detector simulation
class Detector {
  public:
    Detector(runCard const& run, sampleCard const& sample, int const& pythiaSeed,
             int const& detectorSeed);
    void Simulate(finalState in, finalState& out);

  private:
    void AddPileup(finalState& fs);
    Pythia8::Pythia pythiaPileup;
    pcg32           rng;

    const int    nPileup;
    const double phiRes;
    const double etaRes;

    const double jetEsmear;
};