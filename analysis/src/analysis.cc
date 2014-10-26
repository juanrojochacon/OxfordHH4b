// analyis.cc jr/nh 21/10/14

#include "analysis.h"
#include "yodagrams.h"

#include <iostream>

#include "YODA/Histo1D.h"
#include "YODA/WriterFLAT.h"

#include "fastjet/Selector.hh"

  // General string hasher
static int IntHash(const std::string& _str)
{
	const char* s = _str.c_str();
	unsigned h = 31;
	while (*s) {
		h = (h * 54059) ^ (s[0] * 76963);
		s++;
	}
	return h % 86969;
};

Analysis::Analysis(string const& name, string const& sample):
analysisName(name),
analysisRoot("/" + std::string(RESDIR) +"/"+ name + "/"),
sampleName(sample)
{
	std::cout << "Analysis " << analysisName << " initialised at: " <<analysisRoot<<std::endl; 
	outputNTuple.open( "." + analysisRoot + "ntuples/" +sampleName+ "_NTuple.dat");
};


Analysis::~Analysis()
{
	// Close NTuple output
	outputNTuple.close();

	// Export files
	Export();

}


void Analysis::BookHistogram(YODA::Histo1D* hist, string const& name)
{
	// Histo path
	const string path  = analysisRoot + "histodat/"+sampleName+"_" + name + ".dat";

	// Add to histogram prototypes
	hist->setTitle(name);
	hist->setPath(path);

	std::map<int,YODA::Histo1D*>::iterator iMap = bookedHistograms.find(IntHash(name));
	if (iMap != bookedHistograms.end())
	{
		std::cerr << "Analysis::BookHistogram error: HASH COLLISION for histogram: "<<name<<std::endl;
		std::cerr << "Either histogram is duplicated, or we need a better hash function!"<<std::endl;
	}
	else
	{
		bookedHistograms.insert(std::make_pair(IntHash(name),hist));
	}

//	hist->setAnnotation(std::string("XLabel"), std::string("p_T (GeV)"));

}

void Analysis::FillHistogram(string const& rname, double const& weight, double const& coord )
{
	std::map<int,YODA::Histo1D*>::iterator iMap = bookedHistograms.find(IntHash(rname));
	if (iMap != bookedHistograms.end())
	{	(*iMap).second->fill(coord,weight);	}
	else
	{
		std::cerr << "Analysis::FillHistogram error: Cannot find Histogram: "<<rname<<std::endl;
		exit(-1);
	}
}

void Analysis::Export()
{
	// Write out cut flow
	std::cout << "Exporting cutFlow: "<<analysisRoot + "cutFlow.dat"<<std::endl;
	std::ofstream cutFlow("./" + analysisRoot + sampleName + "_cutFlow.dat");
	for (size_t i=0; i<cutWeight.size(); i++)
		cutFlow << i << "\t" << cutWeight[i].first <<"\t" << cutWeight[i].second<<std::endl;

	// Export histograms
	std::map<int,YODA::Histo1D*>::iterator iMap = bookedHistograms.begin();
	while (iMap != bookedHistograms.end())
	{
		std::cout << "Writing Histogram: "<< (*iMap).second->path()<<std::endl;
		YODA::WriterFLAT::write("." + (*iMap).second->path(), *(*iMap).second);
		iMap++;
	}
}

void Analysis::Cut(std::string const& cutStr, double const& weight)
{
	// Check for already booked cuts
	for (size_t i=0; i< cutWeight.size(); i++)
		if (cutWeight[i].first.compare(cutStr) == 0)
		{
			cutWeight[i].second += weight;
			return;
		}

	// No cut found, book a new one
	cutWeight.push_back(std::make_pair<const std::string, double>(cutStr, weight));
	return;
}

void Analysis::Pass(double const& weight)
{
	nPassed++;
	passedWeight += weight;
}