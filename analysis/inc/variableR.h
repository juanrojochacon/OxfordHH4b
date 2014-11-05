// basic.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>

/*
	Basic template analysis
 */
class VariableRAnalysis : public Analysis
{
	public:
		// Constructor for analysis -> books histograms, specified cutFlow etc
		VariableRAnalysis(std::string const& sampleName);

		// Analysis function - takes a finalstate and processes it
		void Analyse(bool const& signal, double const& weightnorm, finalState const&) = 0;

	private:
};