// oxford_combined.h
#pragma once

#include "analysis.h"
#include "utils.h"

#include <iostream>
#include <vector>

/*
combined analysis for HH in fully merged, semi-merged and resolved channels
 */
class OxfordCombinedCheckAnalysis : public Analysis
{
	public:
		OxfordCombinedCheckAnalysis(std::string const& sampleName, int const& subsample);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);
	private:
		// Small-R B-tagging
		void BTagging( std::vector<fastjet::PseudoJet> const& jets_vec, 
					   std::vector<btagType>& btag_vec  );

		// Large-R B-tagging
 		void BTagging( 	std::vector<fastjet::PseudoJet> const& largeRJets, 
 						std::vector<fastjet::PseudoJet> const& trackjets, 
 						std::vector<fastjet::PseudoJet>& subjet1,
 						std::vector<fastjet::PseudoJet>& subjet2,
 						std::vector<int>& nBSubJets_vec,
                      	std::vector<int>& nCSubJets_vec,
                      	std::vector<int>& nLSubJets_vec  );

        // Small-R reconstruction
		void Reco_Resolved( std::vector<fastjet::PseudoJet> const& bjets_vec, 
							std::vector<fastjet::PseudoJet>& higgs_vec,
							std::vector<fastjet::PseudoJet>& higgs0_vec,  	// Leading higgs subjets
                            std::vector<fastjet::PseudoJet>& higgs1_vec  ); // Subleading higgs subjets

		bool Reco_Intermediate(   std::vector<fastjet::PseudoJet> const& bjets, 
	                              fastjet::PseudoJet const& fatjet, 
	                              std::vector<btagType> const& btag_vec,
	                              std::vector<btagType>&  btag_selected_vec,
	                              fastjet::PseudoJet& lead_subjet,
	                              fastjet::PseudoJet& sublead_subjet,
	                              std::vector<fastjet::PseudoJet>& higgs_vec );
};
