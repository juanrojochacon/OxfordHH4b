// basic.cc

#include "basic.h"
#include "utils.h"

// Some useful fastjet includes
//#include "fastjet/Selector.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/tools/MassDropTagger.hh"

#include "YODA/Histo1D.h"

BasicAnalysis::BasicAnalysis(runCard const& run, sampleCard const& sample, int const& subsample)
    : Analysis("basic", run, sample, subsample) {
    // Example histo parameters
    const int    nbins_example = 20;
    const double histo_min     = 0;
    const double histo_max     = 500;

    // Example histogram
    BookHistogram(new YODA::Histo1D(nbins_example, histo_min, histo_max), "example_histo");
}

void BasicAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs) {
    // Call basic analysis
    Analysis::Analyse(signal, weightnorm, fs);

    // Weightnorm provides the sample's unit weight
    double event_weight = weightnorm;

    std::cout << "FinalState: " << std::endl;
    for (size_t i = 0; i < fs.size(); i++)
        std::cout << i << "\t" << fs[i].user_index() << std::endl;
    std::cout << std::endl;

    // Do stuff to the final state
    // DoStuff(fs);

    // Fail a cut
    // if(!passCuts1(fs)) return Cut("cut1", event_weight);
    // Fail a different cut
    // if(!passCuts2(fs)) return Cut("cut2", event_weight);

    // Fill example histogram (after cuts)
    const double hist_coord = 5; // Observable to be binned
    FillHistogram("example_histo", event_weight, hist_coord);

    // Pass event
    Pass(event_weight);
}