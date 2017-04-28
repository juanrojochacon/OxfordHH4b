// basic.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>

/*
    Basic template analysis
 */
class BasicAnalysis : public Analysis {
  public:
    // Constructor for analysis -> books histograms, specified cutFlow etc
    BasicAnalysis(runCard const& run, sampleCard const& sample, int const& subSample);

    // Analysis function - takes a finalstate and processes it
    void Analyse(bool const& signal, double const& weightnorm, finalState const&);

  private:
};