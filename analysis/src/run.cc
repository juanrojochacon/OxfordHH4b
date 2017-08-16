#include "run.h"

#include "analysis.h"
#include "utils.h"

#include "oxford.h"
#include "oxford_atlas_qcd.h"
#include "oxford_sideband.h"
#include "samples.h"

#include <cmath>
#include <sstream>
#include <vector>

// Initialise list of analyses
void InitAnalyses(std::vector<Analysis*>& analyses, runCard const& run, sampleCard const& sample,
                  int const& subsample) {
    // **************** PLEASE MODIFY *****************
    // analyses.push_back(new OxfordAnalysis(run, sample, subsample));
    //  analyses.push_back(new OxfordSidebandAnalysis(run, sample, subsample));
    analyses.push_back(new OxfordAtlasQcdAnalysis(run, sample, subsample, 2));
    analyses.push_back(new OxfordAtlasQcdAnalysis(run, sample, subsample, 3));
    analyses.push_back(new OxfordAtlasQcdAnalysis(run, sample, subsample, 4));
    // **************** DO NOT MODIFY  ****************
}

// ************************* Sample card parsing ********************************

template <class T> T cardquery(std::string const& filename, std::string const& field) {
    std::string   delim = "=";
    std::ifstream instr(filename);
    std::string   line;
    while (getline(instr, line)) {
        const size_t pos = line.find(delim);
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            line.erase(0, pos + delim.length());
            std::string value = line.substr(0, line.length());
            if (key == field) {
                std::stringstream ss;
                T                 outval;
                ss << value;
                ss >> outval;

                instr.close();
                return outval;
            }
        }
    }
    instr.close();
    throw std::runtime_error("Cannot find key: " + field);
}

sampleCard::sampleCard(std::string const& filename)
    : eventfile(cardquery<std::string>(filename, "eventfile")),
      eventpath(std::string(SAMPLEDIR) + eventfile),
      samplename(cardquery<std::string>(filename, "samplename")),
      format(cardquery<std::string>(filename, "format")),
      hepmc(format == "HEPMC"),
      xsec_norm(cardquery<double>(filename, "xsec_norm")),
      sqrts(cardquery<double>(filename, "sqrts")),
      is_signal(cardquery<bool>(filename, "signal")),
      nevt_sample(cardquery<double>(filename, "nevt_sample")) {
    std::cout << "-- Parsed sampleCard ------------------" << std::endl;
    std::cout << "   Sample name: " << samplename << std::endl;
    std::cout << "   Event file:  " << eventfile << std::endl;
    std::cout << "   Format:      " << format << std::endl;
    std::cout << "   xSec norm:   " << xsec_norm << std::endl;
    std::cout << "   CoM energy:  " << sqrts << std::endl;
    std::cout << "   Signal:      " << is_signal << std::endl;
    std::cout << "   N_evt:       " << nevt_sample << std::endl;
    std::cout << "---------------------------------------" << std::endl;
}

runCard::runCard(std::string const& filename)
    : runname(cardquery<std::string>(filename, "runname")),
      sub_samplesize(cardquery<double>(filename, "sub_samplesize")),
      npileup(cardquery<double>(filename, "npileup")),
      jetEsmear(cardquery<double>(filename, "jetEsmear")),
      pythiaShower(cardquery<bool>(filename, "pythiaShower")),
      runseed(cardquery<double>(filename, "runseed")) {
    std::cout << "-- Parsed runCard --------------------" << std::endl;
    std::cout << "   Run name:       " << runname << std::endl;
    std::cout << "   SubSample size: " << sub_samplesize << std::endl;
    std::cout << "   N_PU:           " << npileup << std::endl;
    std::cout << "   Jet E smear:    " << jetEsmear << std::endl;
    std::cout << "   pythiaShower:   " << pythiaShower << std::endl;
    std::cout << "   runseed:        " << runseed << std::endl;
    std::cout << "---------------------------------------" << std::endl;
}
