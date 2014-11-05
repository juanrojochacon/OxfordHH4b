// basic.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>

/*
	Basic template analysis
 */
class BasicAnalysis : public Analysis
{
	public:
		BasicAnalysis(std::string const& sampleName);

		void Analyse(bool const& signal, double const& weightnorm, finalState const&) = 0;

	private:
};