// oxford_atlas_qcd.h
#pragma once

#include "analysis.h"
#include "utils.h"

#include <iostream>
#include <vector>

/*
combined analysis for HH in fully merged, semi-merged and resolved channels - reweighting events
rather than rejecting them at b-tagging
 */
class OxfordAtlasQcdAnalysis : public Analysis {
  public:
    OxfordAtlasQcdAnalysis(runCard const& run, sampleCard const& sample, int const& subSample,
                           int const& nBTag);

    void Analyse(bool const& signal, double const& weight_norm, finalState const&,
                 double gen_weight);

  private:
    double ResolvedAnalysis(std::vector<fastjet::PseudoJet> const& smallRJets,
                            std::vector<btagType> const& btags, bool const& signal,
                            double const& event_weight, double gen_weight);
    double BoostedAnalysis(std::vector<fastjet::PseudoJet> const&              largeRJets,
                           std::vector<std::vector<fastjet::PseudoJet>> const& largeRsubJets,
                           std::vector<std::vector<btagType>> btags, bool const& signal,
                           double const& event_weight, double gen_weight);
    double IntermediateAnalysis(std::vector<fastjet::PseudoJet> const&              largeRJets,
                                std::vector<fastjet::PseudoJet> const&              smallRJets,
                                std::vector<std::vector<fastjet::PseudoJet>> const& largeRsubJets,
                                std::vector<std::vector<btagType>> const&           largeRbtags,
                                std::vector<btagType> const& smallRbtags, bool const& signal,
                                double const& event_weight, double gen_weight);
    // Fill common reconstructed higgs quantities
    void HiggsFill(fastjet::PseudoJet const& H0, fastjet::PseudoJet const& H1,
                   std::string const& analysis, std::string const& ntagcat, size_t const& cut,
                   double const& weight);

    void BoostFill(fastjet::PseudoJet const& H, std::string const& analysis,
                   std::string const& ntagcat, size_t const& cut, double const& weight);

    void BoostFill(fastjet::PseudoJet const& H0, fastjet::PseudoJet const& H1,
                   std::string const& analysis, std::string const& ntagcat, size_t const& cut,
                   double const& weight);

    const bool subtractPU;

    int         m_nBTag;
    std::string m_btag_string;

    std::ofstream resNTuple;
    std::ofstream intNTuple;
    std::ofstream bstNTuple;
    std::ofstream weightsNTuple;
    std::ofstream fullNTuple;
};
