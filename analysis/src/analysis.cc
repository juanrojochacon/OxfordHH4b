// analyis.cc jr/nh 21/10/14

#include "analysis.h"
#include "run.h"

#include <iostream>
#include <iomanip>
#include <sys/stat.h>

#include "YODA/Histo1D.h"
#include "YODA/WriterFLAT.h"

#include "fastjet/Selector.hh"

// Default Verbosity
bool Analysis::Verbose = false;

// Create directory structure
static inline void createPath(std::string path)
{
	mkdir(path.c_str(),0777);
}

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
totalWeight(0),
passedWeight(0)
{
	if (Verbose) std::cout << "Creating Path: "<<"."+analysisRoot<<std::endl;
	createPath("." + analysisRoot);
	createPath("." + analysisRoot + sampleName);

	if (Verbose) std::cout << "Analysis " << analysisName << " initialised at: " <<analysisRoot<<std::endl;
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
		exit(-1);
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
	if (Verbose) std::cout << "Exporting cutFlow: "<<analysisRoot + sampleName + "/cutFlow.dat"<<std::endl;
	const string cutFlow_out = "./" + analysisRoot + sampleName + "/cutFlow.dat";
	std::ofstream cutFlow(cutFlow_out.c_str());

	// Total weights
	double sumWeight = 0.0;
	for (size_t i=0; i<cutWeight.size(); i++)
		sumWeight += cutWeight[i].second;

	cutFlow << std::left<<std::setw(5) << 0 << std::left<<std::setw(25) 
		<< "Total" <<std::left<<std::setw(25) 
		<< totalWeight <<std::left<<std::setw(20)<<std::endl;

	for (size_t i=0; i<cutWeight.size(); i++)
	{
		cutFlow << std::left<<std::setw(5) << i+1 << std::left<<std::setw(25) 
		<< cutWeight[i].first <<std::left<<std::setw(25) 
		<< totalWeight - cutWeight[i].second<<std::left<<std::setw(20)<<std::endl;
		totalWeight -= cutWeight[i].second;
	}
	cutFlow << std::left<<std::setw(5) << cutWeight.size()<<std::left
	<<std::setw(25)<<"(Passed)"<<std::left<<std::setw(25)
	<<passedWeight<<"  ("<<totalWeight<<")"<<std::endl;


	cutFlow.close();

	// Export histograms
	std::map<int,YODA::Histo1D*>::iterator iMap = bookedHistograms.begin();
	while (iMap != bookedHistograms.end())
	{
		if (Verbose) std::cout << "Writing Histogram: "<< (*iMap).second->path()<<std::endl;
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
	cutWeight.push_back(std::make_pair<std::string, double>(cutStr, weight));
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


bool Analysis::VerifyFourMomentum(std::vector<fastjet::PseudoJet> const& jets)
{
	// Smearing breaks four-mom verification
	if (GetPTSmear() > 1E-8) return true;
	// No beam-remnants
	if (!pythiaShowered()) return true; 


	// To check energy-momentum conservation
	double const Eref=14000; // Samples generated for LHC 14 TeV
	double const tol_emom=1.0;

	// Check again four-momentum conservation, this time applied to jets
	// formed from the clustering of quarks and gluons (and beam remnants as well)
	double px_tot=0;
	double py_tot=0;
	double pz_tot=0;
	double E_tot=0;
	for(size_t ij=0;ij<jets.size();ij++)
	{
		px_tot+= jets[ij].px();
		py_tot+= jets[ij].py();
		pz_tot+= jets[ij].pz();
		E_tot+= jets[ij].E();
	}

	// Check energy-momentum conservation
	if( fabs(px_tot) > tol_emom || fabs(py_tot)  > tol_emom 
	|| fabs(pz_tot)  > tol_emom || fabs(E_tot-Eref)  > tol_emom )
	{
		std::cout<<"\n ********************************************************************** \n"<<std::endl;
		std::cout<<"No conservation of energy after jet reconstruction "<<std::endl;
		std::cout<<"N_particles = " << jets.size()<<std::endl;
		std::cout<<"px_tot = "<<px_tot<<std::endl;
		std::cout<<"py_tot = "<<py_tot<<std::endl;
		std::cout<<"pz_tot = "<<pz_tot<<std::endl;
		std::cout<<"E_tot, Eref = "<<E_tot<<" "<<Eref<<std::endl;
		std::cout<<"\n ********************************************************************** \n"<<std::endl;
		return false;
	}
	return true;
}









