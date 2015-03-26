// oxford_comb.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
Showerable Truth-level analysis for HH
 */
class OxfordCombAnalysis : public Analysis
{
	public:
		OxfordCombAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);
	private:

		void AnalyseResolved(bool const& signal, double const& weight_norm, finalState const&);
		void AnalyseBoosted(bool const& signal, double const& weight_norm, finalState const&);

		int BTagging( fastjet::PseudoJet const& jet ) const;
		int CTagging( fastjet::PseudoJet const& jet ) const;
};