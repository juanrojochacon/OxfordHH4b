// analyis.cc jr/nh 21/10/14

#include "analysis.h"
#include "run.h"

#include <sys/stat.h>
#include <functional>
#include <iomanip>
#include <iostream>

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

#include "YODA/WriterFLAT.h"
#include "YODA/WriterYODA.h"

#include "fastjet/Selector.hh"

// Default Verbosity
bool Analysis::Verbose = true;

// Create directory structure
static inline void createPath(const std::string& path) { mkdir(path.c_str(), 0777); }

// General string hasher
static int IntHash(const std::string& _str) {
    std::hash<std::string> str_hash;
    return str_hash(_str);
};

Analysis::Analysis(string const& name, runCard const& run, sampleCard const& sample,
                   int const& subsample)
    : analysisName(name),
      analysisRoot("/" + std::string(RESDIR) + "/" + run.runname + "_" + name + "/"),
      sampleName(sample.samplename),
      subSample(subsample),
      nPassed(0),
      totalWeight(0),
      passedWeight(0) {
    if (Verbose)
        std::cout << "Creating Path: "
                  << "." + analysisRoot << std::endl;
    createPath("." + analysisRoot);
    createPath("." + analysisRoot + sampleName);

    if (Verbose)
        std::cout << "Analysis " << analysisName << " initialised at: " << analysisRoot
                  << std::endl;
};

Analysis::~Analysis() {
    // Export files
    Export();

    // Clear all histograms
    for (auto& iT : bookedHistograms_1D) delete iT.second;
    for (auto& iT : bookedHistograms_2D) delete iT.second;

    bookedHistograms_1D.clear();
    bookedHistograms_2D.clear();
}

void Analysis::BookHistogram(YODA::Histo1D* hist, string const& name) {
    // Add to histogram prototypes
    hist->setTitle(name);
    hist->setPath(analysisRoot + "histo_" + name);

    auto iMap = bookedHistograms_1D.find(IntHash(name));
    if (iMap != bookedHistograms_1D.end()) {
        std::cerr << "Analysis::BookHistogram error: HASH COLLISION for histogram: " << name
                  << std::endl;
        std::cerr << "Attempted to insert: " << name << ", hash: " << IntHash(name) << std::endl;
        std::cerr << "Collided with: " << (*iMap).second->title()
                  << ", hash: " << IntHash((*iMap).second->title()) << std::endl;
        exit(-1);
    }
    else {
        bookedHistograms_1D.insert(std::pair<int, YODA::Histo1D*>(IntHash(name), hist));
    }

    //	hist->setAnnotation(std::string("XLabel"), std::string("p_T (GeV)"));
}

void Analysis::BookHistogram(YODA::Histo2D* hist, string const& name) {
    // Add to histogram prototypes
    hist->setTitle(name);
    hist->setPath(analysisRoot + "histo_" + name);

    auto iMap = bookedHistograms_2D.find(IntHash(name));
    if (iMap != bookedHistograms_2D.end()) {
        std::cerr << "Analysis::BookHistogram error: HASH COLLISION for histogram: " << name
                  << std::endl;
        std::cerr << "Attempted to insert: " << name << ", hash: " << IntHash(name) << std::endl;
        std::cerr << "Collided with: " << (*iMap).second->title()
                  << ", hash: " << IntHash((*iMap).second->title()) << std::endl;
        exit(-1);
    }
    else {
        bookedHistograms_2D.insert(std::pair<int, YODA::Histo2D*>(IntHash(name), hist));
    }

    //	hist->setAnnotation(std::string("XLabel"), std::string("p_T (GeV)"));
}

void Analysis::FillHistogram(string const& rname, double const& weight, double const& coord) {
    auto iMap = bookedHistograms_1D.find(IntHash(rname));
    if (iMap != bookedHistograms_1D.end()) {
        (*iMap).second->fill(coord, weight);
    }
    else {
        std::cerr << "Analysis::FillHistogram error: Cannot find Histogram: " << rname << std::endl;
        exit(-1);
    }
}

void Analysis::FillHistogram(string const& rname, double const& weight, double const& coord1,
                             double const& coord2) {
    auto iMap = bookedHistograms_2D.find(IntHash(rname));
    if (iMap != bookedHistograms_2D.end()) {
        (*iMap).second->fill(coord1, coord2, weight);
    }
    else {
        std::cerr << "Analysis::FillHistogram error: Cannot find Histogram: " << rname << std::endl;
        exit(-1);
    }
}

void Analysis::Export() {
    // Export histograms
    auto iMap1D = bookedHistograms_1D.begin();
    while (iMap1D != bookedHistograms_1D.end()) {
        if (Verbose) std::cout << "Writing Histogram: " << (*iMap1D).second->path() << std::endl;
        if ((*iMap1D).second->numEntries() > 0) {
            std::stringstream path;
            path << analysisRoot + sampleName + "/histo_" + (*iMap1D).second->title() << "."
                 << subSample;
            YODA::WriterYODA::write("." + path.str() + ".yoda", *(*iMap1D).second);
        }
        iMap1D++;
    }

    // Export histograms
    auto iMap2D = bookedHistograms_2D.begin();
    while (iMap2D != bookedHistograms_2D.end()) {
        if (Verbose) std::cout << "Writing Histogram: " << (*iMap2D).second->path() << std::endl;
        if ((*iMap2D).second->numEntries() > 0) {
            std::stringstream path;
            path << analysisRoot + sampleName + "/histo_" + (*iMap2D).second->title() << "."
                 << subSample;
            YODA::WriterYODA::write("." + path.str() + ".yoda", *(*iMap2D).second);
        }
        iMap2D++;
    }
}
