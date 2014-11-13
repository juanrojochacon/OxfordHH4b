// oxford_boost_fr.h
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
class OxfordBoostFRAnalysis : public Analysis
{
	public:
		OxfordBoostFRAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);

	private:
		// Cluster Jets using fixed-R jets
		void JetCluster_LargeFR(fastjet::ClusterSequence const& cs_akt, std::vector<fastjet::PseudoJet>& fatjets, std::vector<double>& split12_vec, std::vector<double>& tau21_vec, double& event_weight);
		
		// Tag bs according to UCL strategy
		bool BTagging(fastjet::PseudoJet const& jet) const; //!< b Tagging method
		bool TwoBTagging(fastjet::PseudoJet const& jet) const; //!< 2b Tagging method
		
};
