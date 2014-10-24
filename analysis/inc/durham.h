// ucl.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
This analysis is inspired by the Durham group, based on jet substructure
See for example the slides of their talk at Boost 2014
https://indico.cern.ch/event/302395/session/12/contribution/25/material/slides/0.pdf
as well as more details in their publication
 arXiv:1404.7139
Here we only use BDRS tagging, no attemp to use Shower Deconstruction
 */
class DurhamAnalysis : public Analysis
{
	public:
		DurhamAnalysis();

		void Analyse(string const& sampleID, bool const& signal, finalState const&);

	private:
		void JetCluster_Durham(finalState const& particles, string const& sampleID, std::vector<fastjet::PseudoJet>& higgs_candidates, double& event_weight);
};