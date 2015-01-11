// basic.cc

#include "basic.h"
#include "utils.h"
#include "settings.h"

// Some useful fastjet includes
//#include "fastjet/Selector.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/tools/MassDropTagger.hh"

#include "YODA/Histo1D.h"
	
BasicAnalysis::BasicAnalysis(string const& sample):
Analysis("basic", sample)
{
	// Example histo parameters
	const int nbins_example = 20;
	const double histo_min = 0;
	const double histo_max = 500;

	// Example histogram
	BookHistogram(new YODA::Histo1D(nbins_example, histo_min, histo_max), "example_histo");

	// Output nTuple kinematics description (# signal source) is required
	const std::string tupleSpec = "# signal source example_kin example_kin2";
 	outputNTuple<<tupleSpec<<std::endl;	// DO NOT CHANGE
}

void BasicAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
	// Call basic analysis
	Analysis::Analyse(signal, weightnorm, fs);

	// Weightnorm provides the sample's unit weight
	double event_weight = weightnorm;

	// Do stuff to the final state
	// DoStuff(fs);

	// Fail a cut
	//if(!passCuts1(fs)) return Cut("cut1", event_weight);
	// Fail a different cut
	//if(!passCuts2(fs)) return Cut("cut2", event_weight);

	// Fill example histogram (after cuts)
	const double hist_coord = 5; // Observable to be binned
	FillHistogram("example_histo", event_weight, hist_coord );

	
	// Write kinematics to nTuple
	outputNTuple <<signal <<"\t"<< GetSample() <<"\t" // DO NOT CHANGE
	<< hist_coord<<"\t"<<5*hist_coord<<"\t"<<std::endl; 

	// Pass event
	Pass(event_weight);
}