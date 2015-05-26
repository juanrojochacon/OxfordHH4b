// oxford_combined.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
combined analysis for HH in fully merged, semi-merged and resolved channels - reweighting events rather than rejecting them at b-tagging
 */
class OxfordCombinedRWAnalysis : public Analysis
{
	public:
		OxfordCombinedRWAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);
	private:
		void BTagging( std::vector<fastjet::PseudoJet> const& jets_vec, std::vector<bool>& isFake_vec  );
 		void BTaggingFJ( std::vector<fastjet::PseudoJet> const& largeRJets, std::vector<fastjet::PseudoJet> const& trackjets, std::vector<int>& nBSubJets_vec );


		void Reco_Resolved( std::vector<fastjet::PseudoJet> const& bjets_vec, std::vector<fastjet::PseudoJet>& higgs_vec );
		bool Reco_Intermediate( std::vector<fastjet::PseudoJet> const& bjets, std::vector<bool> const& isFakeSR_vec, fastjet::PseudoJet const& fatjet, int& nBjets, std::vector<fastjet::PseudoJet>& higgs_vec );

		std::ofstream resNTuple;
		std::ofstream intNTuple;
		std::ofstream bstNTuple;
};
