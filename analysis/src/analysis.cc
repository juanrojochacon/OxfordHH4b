// analyis.cc jr/nh 21/10/14

#include "analysis.h"

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
sampleName(sample),
nPassed(0),
passedWeight(0)
{
	std::cout << "Analysis " << analysisName << " initialised at: " <<analysisRoot<<std::endl;
	const string ntupOut =  "." + analysisRoot + sampleName + "/ntuple.dat";
	outputNTuple.open( ntupOut.c_str() );
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
	const string path  = analysisRoot +sampleName+"/histo_" + name + ".dat";

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
	std::cout << "Exporting cutFlow: "<<analysisRoot + sampleName + "/cutFlow.dat"<<std::endl;
	const string cutFlow_out = "./" + analysisRoot + sampleName + "/cutFlow.dat";
	std::ofstream cutFlow(cutFlow_out.c_str());

	// Total weights
	double sumWeight = 0.0;
	for (size_t i=0; i<cutWeight.size(); i++)
		sumWeight += cutWeight[i].second;

	for (size_t i=0; i<cutWeight.size(); i++)
		cutFlow << i << "\t" << cutWeight[i].first <<"\t" << cutWeight[i].second<<"\t"<< 100.0*(cutWeight[i].second/sumWeight) <<std::endl;
	cutFlow << cutWeight.size()<<"\t"<<"(Passed)"<<"\t"<<passedWeight<<std::endl;
	cutFlow.close();

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

double Analysis::GetCutWeight() const
{
	double cutWeights = 0.0;
	for (size_t i=0; i< cutWeight.size(); i++)
		cutWeights += cutWeight[i].second;

	return cutWeights;
}

void Analysis::Pass(double const& weight)
{
	nPassed++;
	passedWeight += weight;
}