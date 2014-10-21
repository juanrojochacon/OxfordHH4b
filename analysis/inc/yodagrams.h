// yodagrams.h 17/10/14 nh
#pragma once

#include <iostream>

namespace YODA {
	class Histo1D;
}

// Export yoda histograms
void yoda_export();

// Add yoda histogram
void yoda_add(YODA::Histo1D* hist);

// Fill yoda histogram
void yoda_fill(std::string const& histofill, double const& weight, double const& coord);
