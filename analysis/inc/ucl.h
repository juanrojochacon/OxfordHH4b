// ucl.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
This is the analysis used by the UCL group
See for example the slides of their talk at Boost 2014
https://indico.cern.ch/event/302395/session/12/contribution/26/material/slides/1.pdf
 */
class UCLAnalysis : public Analysis
{
	public:
		UCLAnalysis();

		void Analyse(string const& sampleID, bool const& signal, finalState const&);

	private:
		// Cluster Jets according to UCL strategy
		void JetCluster_UCL(finalState const& particles, std::vector<fastjet::PseudoJet>& bjets, double& event_weight);
		
		// Tag bs according to UCL strategy
		bool BTagging(fastjet::PseudoJet const& jet) const; //!< b Tagging method 

};