// oxford_res_fr.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
	This is a resolved style analysis intended to test the b-Tagging procedure
 */
class bTagTestAnalysis : public Analysis
{
	public:
		bTagTestAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);

	private:
		// Test bTagging
		bool BTagging(fastjet::PseudoJet const& jet) const; //!< b Tagging method 
};