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

int main() 
{  
  // Results output
  ofstream out_results;
  const string outDir = "./" + std::string(RESDIR) + "/hh4b.res";
  out_results.open(outDir.c_str()); 

  // Scientific output
  out_results << std::scientific;
  cout << std::scientific;

  // Set here the path to the MC samples Dropbox folder
  string samples_path=std::string(SAMPLEDIR);

  // Initialise analyses
  vector<Analysis*> HH4bAnalyses;
  vector<Analysis*> signalAnalyses;
  vector<Analysis*> backgroundAnalyses;

  // Initialise all-sample analyses
  InitAnalyses(HH4bAnalyses, signalAnalyses, backgroundAnalyses);
  for(int is=0; is<GetNSamples(); is++)
  {
    // Read sample data
    eventSample sample = GetSample(is);

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
    InitSampleAnalyses(sampleAnalyses, sample.samplename);

    // total xsec counter
    double sample_xsec = 0;
    // Begin loop over events
    for (int iEvent = 0; iEvent < sample.nevt_sample; ++iEvent) 
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

      // Total analyses
      for (size_t i=0; i<HH4bAnalyses.size(); i++)
        HH4bAnalyses[i]->Analyse(sample.signal, event_weight, fs);

      // Sample analyses
      for (size_t i=0; i<sampleAnalyses.size(); i++)
        sampleAnalyses[i]->Analyse(sample.signal, event_weight, fs);

      // Signal
      if (sample.signal)
        for (size_t i=0; i<signalAnalyses.size(); i++)
          signalAnalyses[i]->Analyse(sample.signal, event_weight, fs);
      else // Background
        for (size_t i=0; i<backgroundAnalyses.size(); i++)
          backgroundAnalyses[i]->Analyse(sample.signal, event_weight, fs);

    // ****************************** Analyses ******************************

    }

    // Compute the total weight of the sample
    // Should coincide with nev_pass for btag prob of 1.0 and light jet mistag prob of 0.0
    for (size_t i=0; i<sampleAnalyses.size(); i++)
    {
      const double wgt_pass = sampleAnalyses[i]->GetWeight();
      const int nev_pass = sampleAnalyses[i]->GetNPassed();
      const double notCounted = sample_xsec - (sampleAnalyses[i]->GetCutWeight() + sampleAnalyses[i]->GetWeight());

      // Save results for cross-sections and number of events
      // Use LHC Run II and HL-LHC luminosities
      out_results<<"\nSample = "<< sample.samplename<<" , Analysis = "<< sampleAnalyses[i]->GetName()<<std::endl;
      out_results<<"nev_tot(MC), nev_pass(MC) = "<<sample.nevt_sample<<" , "<<nev_pass<<std::endl;
      out_results<<"xsec_tot, xsec_pass (fb) = "<<sample_xsec<< " , "<<wgt_pass<<std::endl;
      out_results << "cutXsec: "<< sampleAnalyses[i]->GetCutWeight() << ", Unaccounted for weights: "<< notCounted <<endl;
      out_results<<"pass weight (300 1/fb) = "<< lumi_run2*wgt_pass<<std::endl; // LHC run II numbers
      out_results<<"pass weight, nev_pass (3000 1/fb) = "<< lumi_hllhc*wgt_pass<<std::endl; // HL-LHC numbers
    }

    // Free sample analyses
    for (size_t i=0; i<sampleAnalyses.size(); i++)
      delete sampleAnalyses[i];

    // Close stream
    hepmc_is.close();
  }

// Finish up
out_results.close();
for (size_t i=0; i<HH4bAnalyses.size(); i++)
  delete HH4bAnalyses[i];

for (size_t i=0; i<signalAnalyses.size(); i++)
  delete signalAnalyses[i];

for (size_t i=0; i<backgroundAnalyses.size(); i++)
  delete backgroundAnalyses[i];

// End of the main program
return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
