// amcatnlo.cc

#include "amcatnlo.h"
#include "utils.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <fastjet/Selector.hh>

#include "YODA/Histo1D.h"

AMCAnalysis::AMCAnalysis(std::string const& sampleName, int const& subsample)
    : Analysis("aMC@NLO", sampleName, subsample) {
    BookHistogram(new YODA::Histo1D(6, 0, 6), "cutflow");
}
void AMCAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs) {
    Analysis::Analyse(signal, weightnorm, fs);
    double event_weight = weightnorm;
    FillHistogram("cutflow", event_weight, 0.1);

    const double             jetR = 0.5;
    fastjet::JetDefinition   akt(fastjet::antikt_algorithm, jetR);
    fastjet::ClusterSequence cs_akt(fs, akt);

    const fastjet::Selector               rap_sel = fastjet::SelectorAbsRapMax(4.0);
    const std::vector<fastjet::PseudoJet> jets_akt_80 =
          sorted_by_pt(rap_sel(cs_akt.inclusive_jets(80.0)));
    const std::vector<fastjet::PseudoJet> jets_akt_100 =
          sorted_by_pt(rap_sel(cs_akt.inclusive_jets(100.0)));

    // Need at least 4 jets with pT > 80
    if ((int)jets_akt_80.size() < 4) return Cut(event_weight);

    // Need at least 1 jet with pT > 100
    if ((int)jets_akt_100.size() < 1) return Cut(event_weight);

    // b-Tagging
    int NbJet_80  = 0;
    int NbJet_100 = 0;

    const double eta_amc_cut_bjet = 2.5;
    for (size_t i = 0; i < jets_akt_80.size(); i++)
        if (fabs(jets_akt_80[i].eta()) < eta_amc_cut_bjet)
            if (BTagging(jets_akt_80[i])) NbJet_80++;
    for (size_t i = 0; i < jets_akt_100.size(); i++)
        if (fabs(jets_akt_100[i].eta()) < eta_amc_cut_bjet)
            if (BTagging(jets_akt_100[i])) NbJet_100++;

    // Require 2 b-Jets at 80GeV
    if (NbJet_80 < 2) return Cut(event_weight);

    // Require 1 b-Jet at 100GeV
    if (NbJet_100 < 1) return Cut(event_weight);

    FillHistogram("cutflow", event_weight, 5.1);
    Pass(event_weight);
}

bool AMCAnalysis::BTagging(fastjet::PseudoJet const& jet) const {
    const std::vector<fastjet::PseudoJet>& jet_constituents = jet.constituents();
    for (size_t i = 0; i < jet_constituents.size(); i++) {
        const int userid = jet_constituents.at(i).user_index();
        if (abs(userid) == 5) return true;
    }
    return false;
}
