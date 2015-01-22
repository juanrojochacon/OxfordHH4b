// oxford_boost_fr.h
#pragma once

#include "analysis.h"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>


/*
This is a minimal analysis with only jet acceptance cuts (as of 14th November 2014)
 */
class OxfordBoostFRAnalysis : public Analysis
{
	public:
		OxfordBoostFRAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);

	private:
		// Cluster Jets using VariableR jets
		void JetCluster_LargeFR(finalState const& fs, std::vector<fastjet::PseudoJet>& fatjets, std::vector<double>& split12_vec, std::vector<double>& tau21_vec, double& event_weight);
		
		// Tag bs 
		int BTagging(fastjet::PseudoJet const& jet) const; //!< 2b Tagging method
		int CTagging(fastjet::PseudoJet const& jet) const; //!< 2b Tagging method
		
};
