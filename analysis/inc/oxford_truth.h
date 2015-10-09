// oxford_truth.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
Truth-level analysis for HH
 */
class OxfordTruthAnalysis : public Analysis
{
	public:
		OxfordTruthAnalysis(std::string const& sampleName, int const& subsample);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);
	private:
};