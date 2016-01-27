/////////////////////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>

#include "run.h"
#include "analysis.h"
#include "samples.h"
#include "detector.h"

using namespace Pythia8;

////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) 
{  

  if (argc != 4)
  {
    cerr << "Error: Wrong number of arguments!"<<endl;
    cerr << "Usage: HH4b <run card> <sample card> <subsample>" <<endl;
    exit(-1);
  }

  // Read run card
  const std::string runfile = std::string(argv[1]);
  const runCard run(runfile);

  // Read sample card 
  const std::string samplefile = std::string(argv[2]);
  const sampleCard sample(samplefile);

  // Determine subsample constants
  const int subsample = atoi(argv[3]);
  const int sampleStart = subsample*sampleSize(); // start point of the subsample

  pythiaSeed() = 382;
  systemSeed() = 2093;

  cout << "Processing sample: " <<sample.samplename<< ", subsample: : "<<subsample;
  cout << ". RNG Seeds - Pythia: " << pythiaSeed() <<". System: "<< systemSeed() <<"."<<endl;

  // Read sample data
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
  cout << "Skipping to startpoint: " << sampleStart <<endl;
  double dum; finalState dum2;
  for (int iEvent = 0; iEvent < sampleStart; ++iEvent) 
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
  const int targetSize = min(sampleSize(), sample.nevt_sample - sampleStart);
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

    if (iEvent % 1000 == 0 )
      cout << iEvent <<" events analysed"<<endl;

    // Normalise event weight
    event_weight *= weight_norm;

    // Add to sample xsec
    sample_xsec += event_weight;

  // ****************************** Analyses ******************************

    // Sample analyses
    for (size_t i=0; i<sampleAnalyses.size(); i++)
      sampleAnalyses[i]->Analyse(sample.is_signal, event_weight, fs);

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
