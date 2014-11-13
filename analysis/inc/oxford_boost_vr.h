// oxford_boost_vr.h
#pragma once

#include "analysis.h"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>


/*
This is the analysis used by the UCL group
with AKT5 jets replaced by small VR jets with Rmax=0.5, Rmin=0.1 and rho=40GeV
See for example the slides of their talk at Boost 2014
https://indico.cern.ch/event/302395/session/12/contribution/26/material/slides/1.pdf
 */
class OxfordBoostVRAnalysis : public Analysis
{
	public:
		OxfordBoostVRAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);

	private:
		// Cluster Jets using VariableR jets
		void JetCluster_LargeVR(fastjet::ClusterSequence const& cs, std::vector<fastjet::PseudoJet>& fatjets, double& event_weight);
		
		// Tag bs 
		bool TwoBTagging(fastjet::PseudoJet const& jet) const; //!< 2b Tagging method
		
};