// oxford_combined.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
combined analysis for HH in fully merged, semi-merged and resolved channels - reweighting events rather than rejecting them at b-tagging
 */
class OxfordCombinedRW2Analysis : public Analysis
{
	public:
		OxfordCombinedRW2Analysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);
	private:
		// B-tagging
		void BTagging( 	std::vector<fastjet::PseudoJet> const& jets_vec, std::vector<bool>& isFake_vec  );
 		void BTagging( 	std::vector<fastjet::PseudoJet> const& largeRJets, 
 						std::vector<fastjet::PseudoJet> const& trackjets, 
 						std::vector<int>& nBSubJets_vec );

 		// Small-R reconstruction
		void Reco_Resolved( std::vector<fastjet::PseudoJet> const& bjets_vec, 
							std::vector<fastjet::PseudoJet>& higgs_vec,
							std::vector<fastjet::PseudoJet>& higgs0_vec,  	// Leading higgs subjets
                            std::vector<fastjet::PseudoJet>& higgs1_vec  ); // Subleading higgs subjets

		bool Reco_Intermediate( std::vector<fastjet::PseudoJet> const& bjets, 
								std::vector<bool> const& isFakeSR_vec, 
								fastjet::PseudoJet const& fatjet, int& nBjets, 
								std::vector<fastjet::PseudoJet>& higgs_vec );

		// Fill basic jet quantities
		void JetFill( 	std::vector<fastjet::PseudoJet> const& jets,
						std::string const& analysis, 
	                    size_t const& cut, 
	                    double const& weight );

		// Fill common reconstructed higgs quantities
		void HiggsFill(	fastjet::PseudoJet const& H0,
	                    fastjet::PseudoJet const& H1,
	                    std::string const& analysis, 
	                    size_t const& cut, 
	                    double const& weight);

		void BoostFill( fastjet::PseudoJet const& H,
						std::string const& analysis, 
	                    size_t const& cut, 
	                    double const& weight );

		void BoostFill( fastjet::PseudoJet const& H0,
						fastjet::PseudoJet const& H1,
						std::string const& analysis, 
	                    size_t const& cut, 
	                    double const& weight );	

		std::ofstream resNTuple;
		std::ofstream intNTuple;
		std::ofstream bstNTuple;
};
