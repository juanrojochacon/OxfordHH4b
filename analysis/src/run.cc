#include "run.h"

#include "utils.h"
#include "analysis.h"

#include "oxford.h"
#include "samples.h"

#include <vector>
#include <cmath>
#include <sstream>

// **************** PLEASE MODIFY  ****************

// Global run parameters
const double jetp_smear = 5.0; // % smear on jet momentum
const double jetE_smear = 5.0; // % smear on jet energy

const bool pythiaShower = true; // Shower events in pythia
const bool softKiller = true; // Enable softkiller pileup removal
const int npileup = 80; // Number of pileup events per hard event

const int samplesize = 3E4; //!< Size of individual subsamples

// **************** DO NOT MODIFY  ****************

static double random_seed_pythia = 40487;	//!< Random seed for pythia
static double random_seed_system = 23429;

double GetPSmear() {return jetp_smear;};
double GetESmear() {return jetE_smear;};

bool pythiaShowered() {return pythiaShower;};
bool softKillered() {return softKiller;};

int sampleSize() {return samplesize;};

double& pythiaSeed() {return random_seed_pythia;};
double& systemSeed() {return random_seed_system;};

int npileupEvents() {return npileup;};


void InitSampleAnalyses( std::vector<Analysis*>& sampleAnalyses, std::string const& samplename, int const& subsample )
{
	// **************** PLEASE MODIFY *****************
	sampleAnalyses.push_back(new OxfordAnalysis(samplename, subsample));

	// **************** DO NOT MODIFY  ****************
}

// ************************* Sample card parsing ********************************


template<class T>
T cardquery(std::string const& filename, std::string const& field)
{
  std::string delim = "=";
  std::ifstream instr(filename);
  std::string line;
  while(getline(instr, line))
  {
    const size_t pos =  line.find(delim);
    if (pos != std::string::npos)
    {
      std::string key = line.substr(0, pos); 
      line.erase(0, pos + delim.length());
      std::string value = line.substr(0, line.length());
      if (key == field)
      {
        std::stringstream ss;
        T outval;
        ss << value;
        ss >> outval;

        instr.close();
        return outval;
      }
    }
  }
  instr.close();
  throw std::runtime_error("Cannot find key: " + field);
}

  sampleCard::sampleCard(std::string const& filename):
  eventfile(cardquery<std::string>(filename,"eventfile")),
  samplename(cardquery<std::string>(filename,"samplename")),
  format(cardquery<std::string>(filename,"format")),
  hepmc(format == "HEPMC"),
  xsec_norm(cardquery<double>(filename,"xsec_norm")),
  sqrts(cardquery<double>(filename,"sqrts")),
  is_signal(cardquery<bool>(filename,"signal")),
  nevt_sample(cardquery<double>(filename,"nevt_sample"))
  {
    std::cout << "-- Parsed sampleCard ------------------"<<std::endl;
    std::cout << "   Sample name: "<<samplename <<std::endl;
    std::cout << "   Event file:  "<<eventfile <<std::endl;
    std::cout << "   Format:      "<<format <<std::endl;
    std::cout << "   xSec norm:   "<<xsec_norm <<std::endl;
    std::cout << "   CoM energy:  "<<sqrts <<std::endl;
    std::cout << "   Signal:      "<<is_signal <<std::endl;
    std::cout << "   N_evt:       "<<nevt_sample <<std::endl;
    std::cout << "---------------------------------------"<<std::endl;
  }

  runCard::runCard(std::string const& filename):
  runname(cardquery<std::string>(filename, "runname")),
  sub_samplesize(cardquery<double>(filename,"sub_samplesize")),
  npileup(cardquery<double>(filename,"npileup")),    
  jetEsmear(cardquery<double>(filename,"jetEsmear")),      
  pythiaShower(cardquery<bool>(filename,"pythiaShower")),      
  softKillered(cardquery<bool>(filename,"softKillered")),      
  runseed(cardquery<double>(filename,"runseed"))
  {
    std::cout << "-- Parsed runCard --------------------"<<std::endl;
    std::cout << "   Run name:       "<<runname <<std::endl;
    std::cout << "   SubSample size: "<<sub_samplesize <<std::endl;
    std::cout << "   N_PU:           "<<npileup <<std::endl;
    std::cout << "   Jet E smear:    "<<jetEsmear <<std::endl;
    std::cout << "   pythiaShower:   "<<pythiaShower <<std::endl;
    std::cout << "   softKiller:    "<<softKillered <<std::endl;
    std::cout << "   runseed:       "<<runseed <<std::endl;
    std::cout << "---------------------------------------"<<std::endl;
  }
