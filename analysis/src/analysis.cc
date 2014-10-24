// analyis.cc jr/nh 21/10/14

#include "analysis.h"
#include "yodagrams.h"

#include <iostream>

#include "YODA/Histo1D.h"

#include "fastjet/Selector.hh"


Analysis::Analysis(string const& name):
analysisName(name),
analysisRoot("/" + std::string(RESDIR) +"/"+ name + "/")
{
	std::cout << "Analysis " << analysisName << " initialised at: " <<analysisRoot<<std::endl; 
	totalNTuple.open( "." + analysisRoot + "ntuples/totalNTuple.dat");
};


Analysis::~Analysis()
{
	totalNTuple.close();
	if (sampleNTuple.is_open())
		sampleNTuple.close();
}


void Analysis::BookHistogram(YODA::Histo1D* hist, string const& name)
{
	// Add to histogram prototypes
	hist->setTitle(name);
	protoHistograms.push_back(hist);

	// Generate total sample histogram
	const string totalPath  = analysisRoot + "histodat/Total_" + name + ".dat";
	YODA::Histo1D* totalHisto = new YODA::Histo1D(*hist);
	totalHisto->setPath(totalPath);

	// Add histo to analysis
	yoda_add(totalHisto);
	bookedHistograms.push_back(totalPath);
}

void Analysis::FillHistogram(string const& rname, string const& sname, double const& weight, double const& coord )
{
	// Path for both total and sample histograms
	const string totalPath  = analysisRoot + "histodat/Total_" + rname + ".dat";
	const string samplePath = analysisRoot + "histodat/" + sname + "_" + rname + ".dat";

	// Fill histograms
	yoda_fill(totalPath,  weight, coord);
	yoda_fill(samplePath, weight, coord);
}


void Analysis::InitSample(string const& sampleName)
{
	// Generate histograms
	std::list<YODA::Histo1D*>::iterator iProto;
	for (iProto = protoHistograms.begin(); iProto != protoHistograms.end(); iProto++)
	{
		// Generate new sample histogram
		const string samplePath  = analysisRoot + "histodat/" + sampleName + "_" + (*iProto)->title() + ".dat";
		YODA::Histo1D* sampleHisto = new YODA::Histo1D(*(*iProto));
		sampleHisto->setPath(samplePath);

		// Add to yoda list
		yoda_add(sampleHisto);
	}

	// new Ntuple
	if (sampleNTuple.is_open())
		sampleNTuple.close();
	const string ntup_path = analysisRoot + "ntuples/" + sampleName + "NTuple.dat";
	std::cout << "Writing NTuple to " << ntup_path << std::endl; 
	sampleNTuple.open("." + analysisRoot + "ntuples/" + sampleName + "NTuple.dat");
	sampleNTuple<<tupleSpec<<std::endl;
}

void Analysis::ClearWeights()
{
	// Reset event counters
	nPassed = 0;
	passedWeight = 0;
}