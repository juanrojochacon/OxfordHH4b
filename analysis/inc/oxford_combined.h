// oxford_combined.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
combined analysis for HH in fully merged, semi-merged and resolved channels
 */
class OxfordCombinedAnalysis : public Analysis
{
	public:
		OxfordCombinedAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);
	private:
		void BTagging( std::vector<fastjet::PseudoJet>& jets_vec, std::vector<int>& nBQuarks_vec, std::vector<bool>& isBTagged_vec  );
		void BTaggingFJ( double event_weight, std::vector<fastjet::PseudoJet>& largeRJets, std::vector<fastjet::PseudoJet>& trackjets, 
				 std::vector<int>& nSubJets_vec,  std::vector<int>& nBSubJets_vec,  std::vector<int>& nBTaggedSubJets_vec );
		void Reco_Resolved( std::vector<fastjet::PseudoJet>& bjets_vec, std::vector<fastjet::PseudoJet>& higgs_vec );
		bool Reco_Intermediate( std::vector<fastjet::PseudoJet>& bjets, fastjet::PseudoJet& fatjet, std::vector<fastjet::PseudoJet>& higgs_vec );
};
