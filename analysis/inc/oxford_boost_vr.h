// oxford_boost_vr.h
#pragma once

#include "analysis.h"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>


/*
This is a minimal analysis with only jet acceptance cuts (as of 14th November 2014)
 */
class OxfordBoostVRAnalysis : public Analysis
{
	public:
		OxfordBoostVRAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);

	private:
		// Cluster Jets using VariableR jets
		void JetCluster_LargeVR(finalState const& fs, std::vector<fastjet::PseudoJet>& fatjets, std::vector<double>& split12_vec, std::vector<double>& tau21_vec, double& event_weight);
		
		// Tag bs 
		bool BTagging(fastjet::PseudoJet const& jet) const; //!< 2b Tagging method
		bool CTagging(fastjet::PseudoJet const& jet) const; //!< 2b Tagging method
		
};