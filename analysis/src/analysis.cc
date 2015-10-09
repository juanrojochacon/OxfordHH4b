// analyis.cc jr/nh 21/10/14

#include "analysis.h"
#include "run.h"

#include <iostream>
#include <iomanip>
#include <functional>
#include <sys/stat.h>

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

#include "YODA/WriterFLAT.h"
#include "YODA/WriterYODA.h"

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
	std::hash<std::string> str_hash;
	return str_hash(_str);
};

Analysis::Analysis(string const& name, string const& sample, int const& subsample):
analysisName(name),
analysisRoot("/" + std::string(RESDIR) +"/"+ name + "/"),
sampleName(sample),
subSample(subsample),
nPassed(0),
totalWeight(0),
passedWeight(0)
{
	if (Verbose) std::cout << "Creating Path: "<<"."+analysisRoot<<std::endl;
	createPath("." + analysisRoot);
	createPath("." + analysisRoot + sampleName);

	if (Verbose) std::cout << "Analysis " << analysisName << " initialised at: " <<analysisRoot<<std::endl;
};


Analysis::~Analysis()
{
	// Export files
	Export();
}


void Analysis::BookHistogram(YODA::Histo1D* hist, string const& name)
{
	// Add to histogram prototypes
	hist->setTitle(name);
	hist->setPath(analysisRoot + "histo_" + name);

	std::map<int,YODA::Histo1D*>::iterator iMap = bookedHistograms_1D.find(IntHash(name));
	if (iMap != bookedHistograms_1D.end())
	{
		std::cerr << "Analysis::BookHistogram error: HASH COLLISION for histogram: "<<name<<std::endl;
		std::cerr << "Attempted to insert: " << name <<", hash: "<<IntHash(name)<<std::endl;
		std::cerr << "Collided with: " << (*iMap).second->title() <<", hash: "<<IntHash((*iMap).second->title())<<std::endl;
		exit(-1);
	}
	else
	{
		bookedHistograms_1D.insert(std::pair<int,YODA::Histo1D*>(IntHash(name),hist));
	}

//	hist->setAnnotation(std::string("XLabel"), std::string("p_T (GeV)"));

}

void Analysis::BookHistogram(YODA::Histo2D* hist, string const& name)
{
	// Add to histogram prototypes
	hist->setTitle(name);
	hist->setPath(analysisRoot + "histo_" + name);

	std::map<int,YODA::Histo2D*>::iterator iMap = bookedHistograms_2D.find(IntHash(name));
	if (iMap != bookedHistograms_2D.end())
	{
		std::cerr << "Analysis::BookHistogram error: HASH COLLISION for histogram: "<<name<<std::endl;
		std::cerr << "Attempted to insert: " << name <<", hash: "<<IntHash(name)<<std::endl;
		std::cerr << "Collided with: " << (*iMap).second->title() <<", hash: "<<IntHash((*iMap).second->title())<<std::endl;
		exit(-1);
	}
	else
	{
		bookedHistograms_2D.insert(std::pair<int,YODA::Histo2D*>(IntHash(name),hist));
	}

//	hist->setAnnotation(std::string("XLabel"), std::string("p_T (GeV)"));

}

void Analysis::FillHistogram(string const& rname, double const& weight, double const& coord )
{
	std::map<int,YODA::Histo1D*>::iterator iMap = bookedHistograms_1D.find(IntHash(rname));
	if (iMap != bookedHistograms_1D.end())
	{	(*iMap).second->fill(coord,weight);	}
	else
	{
		std::cerr << "Analysis::FillHistogram error: Cannot find Histogram: "<<rname<<std::endl;
		exit(-1);
	}
}

void Analysis::FillHistogram(string const& rname, double const& weight, double const& coord1, double const& coord2 )
{
	std::map<int,YODA::Histo2D*>::iterator iMap = bookedHistograms_2D.find(IntHash(rname));
	if (iMap != bookedHistograms_2D.end())
	{	(*iMap).second->fill(coord1,coord2,weight);	}
	else
	{
		std::cerr << "Analysis::FillHistogram error: Cannot find Histogram: "<<rname<<std::endl;
		exit(-1);
	}
}

void Analysis::Export()
{
	// Export histograms
	std::map<int,YODA::Histo1D*>::iterator iMap1D = bookedHistograms_1D.begin();
	while (iMap1D != bookedHistograms_1D.end())
	{
		if (Verbose) std::cout << "Writing Histogram: "<< (*iMap1D).second->path()<<std::endl;
		if ((*iMap1D).second->numEntries() > 0)
		{
			std::stringstream path;
			path << analysisRoot +sampleName+"/histo_" + (*iMap1D).second->title() << "." <<subSample;
			YODA::WriterFLAT::write("." + path.str() + ".dat", *(*iMap1D).second);
			YODA::WriterYODA::write("." + path.str() + ".yoda", *(*iMap1D).second);
		}
		iMap1D++;
	}

	// Export histograms
	std::map<int,YODA::Histo2D*>::iterator iMap2D = bookedHistograms_2D.begin();
	while (iMap2D != bookedHistograms_2D.end())
	{
		if (Verbose) std::cout << "Writing Histogram: "<< (*iMap2D).second->path()<<std::endl;
		if ((*iMap2D).second->numEntries() > 0)
		{
			std::stringstream path;
			path << analysisRoot +sampleName+"/histo_" + (*iMap2D).second->title() <<"."<<subSample;
			YODA::WriterFLAT::write("." + path.str() + ".dat", *(*iMap2D).second);
			YODA::WriterYODA::write("." + path.str() + ".yoda", *(*iMap2D).second);
		}
		iMap2D++;
	}
}


bool Analysis::VerifyFourMomentum(std::vector<fastjet::PseudoJet> const& jets)
{
	// Smearing breaks four-mom verification
	if (GetPSmear() > 1E-8 || GetESmear() > 1E-8 ) return true;
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
