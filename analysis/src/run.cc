#include "run.h"

#include "settings.h"
#include "utils.h"
#include "analysis.h"

#include "oxford_combined_rw2.h"
#include "oxford_combined_rw2_check.h"

#include <vector>
#include <cmath>
#include <sstream>

// List of samples
static std::vector<eventSample> samples;

	// **************** PLEASE MODIFY  ****************

// Global run parameters
const int nSamples = 5;
const int max_evt = 1E7;

const double jetp_smear = 5.0; // % smear on jet momentum
const double jetE_smear = 5.0; // % smear on jet energy

const bool pythiaShower = true; // Shower events in pythia
const bool softKiller = true; // Enable softkiller pileup removal

const bool pileup = true;
const int npileup = 80;
const int npileup_total = 3999999;

const int samplesize = 3E5; //!< Size of individual subsamples

const std::string minbias = "PYTHIA_MinBias_14TEV.";

// **************** DO NOT MODIFY  ****************

// Current working subsample
static int subsample = -1;

static double random_seed_pythia = 40487;	//!< Random seed for pythia
static double random_seed_system = 23429;

double GetPSmear() {return jetp_smear;};
double GetESmear() {return jetE_smear;};

int GetNSamples() {return nSamples;};
bool pythiaShowered() {return pythiaShower;};
bool softKillered() {return softKiller;};

int& subSample() {return subsample;};
int sampleSize() {return samplesize;};
int sampleStart() {return samplesize*subsample;};

double& pythiaSeed() {return random_seed_pythia;};
double& systemSeed() {return random_seed_system;};

bool pileupSimulated() {return pileup;};
int npileupEvents() {return npileup;};
int npileupTotal() {return npileup_total;};

std::string minBiasFile() 
{
	std::stringstream filename;
	filename << minbias <<subSample() <<".hepmc";
	return filename.str();
};

eventSample GetSample( int const& isample )
{
	// Populate typical info
	string eventfile;
	string samplename;
	double xsec_norm = 1E3; // pb -> fb
	bool signal;
	bool hepmc;
	int nevt_sample;

	switch (isample)
	{

	// **************** PLEASE MODIFY  ****************
	  case 0:
	  eventfile="HH_sm_eft_1M.lhe";
	  samplename="diHiggs";
	  signal = true;
	  hepmc = false;
	  nevt_sample = 1E6;
	  xsec_norm *= 2.4; // NNLO+NNLL K-factor
	  break;

	  case 1: 
	  eventfile="SHERPA_QCD_2b2j.hepmc";
	  samplename="SHERPA_QCD2b2j";
	  signal = false;
	  hepmc = true;
	  nevt_sample = 3E6;
	  xsec_norm *= 1.3; // NLO K-factor
	  break;

	  case 2: 
	  eventfile="SHERPA_QCD_4b.hepmc";
	  samplename="SHERPA_QCD4b";
	  signal = false;
	  hepmc = true;
	  nevt_sample = 3E6;
  	  xsec_norm *= 1.6; // NLO K-factor
	  break;

	  case 3: 
	  eventfile="SHERPA_QCD_4j.hepmc";
	  samplename="SHERPA_QCD4j";
	  signal = false;
	  hepmc = true;
	  nevt_sample = 3E6;
  	  xsec_norm *= 0.5; // NLO K-factor
	  break;

	  case 4: 
	  eventfile="SHERPA_QCD_ttbar.hepmc";
	  samplename="SHERPA_QCDttbar";
	  signal = false;
	  hepmc = true;
	  nevt_sample = 3E6;
  	  xsec_norm *= 1.4; // NNLO+NNLL K-factor
	  break;


	// **************** DO NOT MODIFY  ****************

	  default:
	  std::cout<<"Invalid Monte Carlo sample, exit"<<std::endl;
	  exit(-10);
	  break;
	}

  eventSample newsamp = 
  {
  	eventfile,
  	samplename,
  	xsec_norm,
  	signal,
  	hepmc,
  	static_cast<int>(fmin(nevt_sample, max_evt))
  };

  return newsamp;
}

void InitSampleAnalyses( std::vector<Analysis*>& sampleAnalyses, std::string const& samplename, int const& subsample )
{
	// **************** PLEASE MODIFY *****************

	sampleAnalyses.push_back(new OxfordCombinedCheckAnalysis(samplename, subsample));
	sampleAnalyses.push_back(new OxfordCombinedRW2Analysis(samplename, subsample));

	// **************** DO NOT MODIFY  ****************
}
