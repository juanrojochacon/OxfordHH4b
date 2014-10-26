// oxford.h
#pragma once

#include "analysis.h"

#include <iostream>
#include <vector>


/*
This is the analysis used by the us
 */
class OxfordAnalysis : public Analysis
{
	public:
		OxfordAnalysis();

		void Analyse(string const& sampleID, bool const& signal, finalState const&);

	private:
};