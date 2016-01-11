/////////////////////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/VariableRPlugin.hh"
#include "Pythia8/Pythia.h"
#include "Pythia8/Info.h"

#include "settings.h"
#include "utils.h"
#include "analysis.h"

#include "pythia.h"
#include "hepmc.h"
#include "detector.h"

#include "run.h"

using namespace Pythia8;

////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) 
{  

  if (argc != 3)
  {
    cerr << "Error: Wrong number of arguments!"<<endl;
    cerr << "Usage: HH4b <sampleID> <subsample>" <<endl;
    exit(-1);
  }

  // Read sampleID
  const int sampleID = atoi(argv[1]);
  const int subsample = atoi(argv[2]);
  int minbiasSample = subsample;
  if( sampleID>0) minbiasSample = subsample+sampleID*100-66;

  subSample() = subsample;
  minBiasSample() = minbiasSample;
  pythiaSeed() = 430598*sampleID +342*subsample + 382;
  systemSeed() = 175*sampleID +34562*subsample + 2093;

  cout << "Processing sample ID: " <<sampleID<< ", subsample: : "<<subsample;
  cout << "MinBias sample ID: " <<minBiasSample();
  cout << ". RNG Seeds - Pythia: " << pythiaSeed() <<". System: "<< systemSeed() <<"."<<endl;

  // Read sample data
  const eventSample sample = GetSample(sampleID);
  const string samples_path=std::string(SAMPLEDIR);
  const string eventfile = samples_path + sample.eventfile;
  std::cout << "Reading samples from: "<<eventfile<<std::endl;

  // Initialise Pythia and HepMC
  Pythia pythiaRun(std::string(PYTHIADIR)); // Pythia input
  std::ifstream hepmc_is( eventfile.c_str() );      // HepMC input

  // Initialise the event sample and weight normalisation
  double weight_norm = 0;
  if (!sample.hepmc)
   InitPythia(pythiaRun, eventfile, sample.nevt_sample, weight_norm );
  else
   InitHepMC( eventfile, sample.nevt_sample, weight_norm );

 // normalise by cross-section normalisation
  weight_norm *= sample.xsec_norm;

  // Initialse Analyses for sample
  vector<Analysis*> sampleAnalyses;
  InitSampleAnalyses(sampleAnalyses, sample.samplename, subsample);

  // Skip to subsample x
  cout << "Skipping to startpoint: " << sampleStart() <<endl;
  double dum; finalState dum2;
  for (int iEvent = 0; iEvent < sampleStart(); ++iEvent) 
  {
    dum2.clear();
    if (!sample.hepmc) // Pythia
      get_final_state_particles(pythiaRun, dum2, dum);
    else  // HepMC      
      get_final_state_particles(hepmc_is,  dum2, dum);
  }

  // total xsec counter
  double sample_xsec = 0;
  // Begin loop over events
  cout << "*************** Analysis Begins ***************" <<endl;
  const int targetSize = min(sampleSize(), sample.nevt_sample - sampleStart());
  cout << "Analysing: " << targetSize <<" events"<<endl;
  for (int iEvent = 0; iEvent < targetSize; ++iEvent) 
  {
    finalState ifs, fs; // The event final state

    double event_weight=0;
    if (!sample.hepmc) // Pythia
      get_final_state_particles(pythiaRun, ifs, event_weight);
    else  // HepMC      
      get_final_state_particles(hepmc_is,  ifs, event_weight);

    // Perform detector simulation
    DetectorSim(ifs,fs);

    if (iEvent % 1000 == 0 && sample.hepmc)
      cout << iEvent <<" HepMC events analysed"<<endl;

    // Normalise event weight
    event_weight *= weight_norm;

    // Add to sample xsec
    sample_xsec += event_weight;

  // ****************************** Analyses ******************************

    // Sample analyses
    for (size_t i=0; i<sampleAnalyses.size(); i++)
      sampleAnalyses[i]->Analyse(sample.signal, event_weight, fs);

  // ****************************** Analyses ******************************

  }

  // Free sample analyses
  for (size_t i=0; i<sampleAnalyses.size(); i++)
    delete sampleAnalyses[i];

  // Close stream and finish up
  hepmc_is.close();

  // End of the main program
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
