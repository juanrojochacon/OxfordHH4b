// oxford.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
Basic analysis
 */
class BasicAnalysis : public Analysis
{
	public:
		BasicAnalysis( string const& sample );

		void Analyse(bool const& signal, double const& weightnorm, finalState const&);

	private:
};