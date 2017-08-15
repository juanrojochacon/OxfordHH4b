// oxford_combined.h
#pragma once

#include "analysis.h"
#include "utils.h"

#include <iostream>
#include <vector>

/*
combined analysis for HH in fully merged, semi-merged and resolved channels - reweighting events
rather than rejecting them at b-tagging
 */
class OxfordAnalysis : public Analysis {
  public:
    OxfordAnalysis(runCard const& run, sampleCard const& sample, int const& subSample);

    void Analyse(bool const& signal, double const& weight_norm, finalState const&, double gen_weight);

  private:
    // Small-R B-tagging
    void BTagging(std::vector<fastjet::PseudoJet> const& jets_vec, std::vector<btagType>& btag_vec);

    // Large-R B-tagging
    void BTagging(std::vector<fastjet::PseudoJet> const& largeRJets,
                  std::vector<fastjet::PseudoJet> const& trackjets,
                  std::vector<fastjet::PseudoJet>&       subjet1,
                  std::vector<fastjet::PseudoJet>& subjet2, std::vector<int>& nBSubJets_vec,
                  std::vector<int>& nCSubJets_vec, std::vector<int>& nLSubJets_vec);

    // Small-R reconstruction
    void Reco_Resolved(std::vector<fastjet::PseudoJet> const& bjets_vec,
                       std::vector<fastjet::PseudoJet>&       higgs_vec,
                       std::vector<fastjet::PseudoJet>&       higgs0_vec, // Leading higgs subjets
                       std::vector<fastjet::PseudoJet>& higgs1_vec); // Subleading higgs subjets

    bool Reco_Intermediate(std::vector<fastjet::PseudoJet> const& bjets,
                           fastjet::PseudoJet const& fatjet, std::vector<btagType> const& btag_vec,
                           std::vector<btagType>& btag_selected_vec,
                           fastjet::PseudoJet& lead_subjet, fastjet::PseudoJet& sublead_subjet,
                           std::vector<fastjet::PseudoJet>& higgs_vec);

    // Fill basic jet quantities
    void JetFill(std::vector<fastjet::PseudoJet> const& smallRJets,
                 std::vector<fastjet::PseudoJet> const& largeRJets, std::string const& analysis,
                 size_t const& cut, double const& weight);

    // Fill basic jet quantities
    void SubJetFill(std::vector<fastjet::PseudoJet> const& leading_subjet,
                    std::vector<fastjet::PseudoJet> const& subleading_subjet,
                    std::string const& analysis, size_t const& cut, double const& weight);

    // Fill common reconstructed higgs quantities
    void HiggsFill(fastjet::PseudoJet const& H0, fastjet::PseudoJet const& H1,
                   std::string const& analysis, size_t const& cut, double const& weight);

    void BoostFill(fastjet::PseudoJet const& H, std::string const& analysis, size_t const& cut,
                   double const& weight);

    void BoostFill(fastjet::PseudoJet const& H0, fastjet::PseudoJet const& H1,
                   std::string const& analysis, size_t const& cut, double const& weight);

    const bool subtractPU;

    std::ofstream resNTuple;
    std::ofstream intNTuple;
    std::ofstream bstNTuple;
};
