// oxford_combined.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
 */
class TestAnalysis : public Analysis
{
	public:
		TestAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weight_norm, finalState const&);
	private:
};
