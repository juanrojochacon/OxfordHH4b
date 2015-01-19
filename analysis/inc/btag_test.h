// oxford_res_fr.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
	This is a resolved style analysis intended to test the b-Tagging procedure
	In this case, only the hardest 4 jets are considered
 */
class bTagTestAnalysis : public Analysis
{
	public:
		bTagTestAnalysis(std::string const& sampleName, std::string analysisName = "bTagTest_hardest4");

		virtual void Analyse(bool const& signal, double const& weight_norm, finalState const&);

	protected:
		// Test bTagging
		int BTagging(fastjet::PseudoJet const& jet) const; //!< b Tagging method - returns number of consitutent b's
};


/*
	This is a resolved style analysis intended to test the b-Tagging procedure
	In this case, all jets with pT > 40 GeV are considered
 */
class bTagTestUCLAnalysis : public bTagTestAnalysis
{
	public:
		bTagTestUCLAnalysis(std::string const& sampleName);

		virtual void Analyse(bool const& signal, double const& weight_norm, finalState const&);
};
