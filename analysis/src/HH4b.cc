/////////////////////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>

// FastJet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
// FastJet contrib
#include "fastjet/contrib/VariableRPlugin.hh"

// Pythia8
#include "Pythia8/Pythia.h"

#include "utils.h"
#include "analysis.h"
#include "ucl.h"
#include "durham.h"

#include "settings.h"

using namespace Pythia8;

////////////////////////////////////////////////////////////////////////////////////////

int main() 
{  

  // Results output
  ofstream out_results;
  out_results.open("./" + std::string(RESDIR) + "/hh4b.res");

  // Set here the path to the MC samples Dropbox folder
  string samples_path=std::string(SAMPLEDIR);

  // Initialise analyses
  vector<Analysis*> HH4bAnalyses;
  HH4bAnalyses.push_back(new UCLAnalysis("total"));
  HH4bAnalyses.push_back(new DurhamAnalysis("total"));
  
  /* ---------------------------------------------------------------------------
  //
  // Loop over signal and background events
  // Signal: gg -> hh with mg5_amc at LO in EFT with heavy quark mass effects
  // Background: QCD 4b, 2b2j, ttbar fully hadronic
  // Need to add new sources of background, in particular QCD 4j
  // Other backgrounds, like single Higgs, also need to be included
  //
  ---------------------------------------------------------------------------*/

  int const nlhe=6; // Number of signal and bkg MC samples

  for(int ilhe=0; ilhe<nlhe;ilhe++){

    // Name of the input file
    string eventfile;
    string samplename;
    bool signal = false;
    
    // Now we loop over the MC samples
    // The total cross-section for each sample can be
    // read from the .lhe file, see the number next to:
    // #  Integrated weight (pb)  
    double xsec = 0.9;

    switch (ilhe)
    {

    // SM gg->HH, 100K
      case 0: 
      eventfile="HH_sm_eft_100K.lhe";
      samplename="SMggHH";
      signal = true;
      xsec = 0.17291e-1;// pb
      xsec *= 1e3; // fb
      xsec *= 0.34; // HH-> 4b BR
      xsec *= 2.5; // NNLO K-factor
      break;

    // BSM gg->HH, 100K, trilinear multiplied by a factor 10
      case 1: 
      eventfile="HH_bsm_lam_10_eft_100K.lhe";
      samplename="BSMggHH";
      signal = true;
      xsec = .27923E+00;// pb
      xsec *= 1e3; // fb
      xsec *= 0.34; // HH-> 4b BR
      xsec *= 2.5; // NNLO K-factor
      break;

    // QCD background, 4b
      case 2: 
      eventfile="qcd_madgraph_4b_14tev_100k_gcuts.lhe";
      samplename="QCD4b";
      signal = false;
      xsec = .58013e+03; // pb
      xsec *= 1e3; // fb
      // Need to add NLO K-factor here
      break;

    // QCD background, 2b2j
      case 3: 
      eventfile="qcd_madgraph_2b2j_14tev_100k_gcuts.lhe";
      samplename="QCD2b2j";
      signal = false;
      xsec = .53822E+06; // pb
      xsec *= 1e3; // fb
      // Need to add NLO K-factor here
      break;

    // QCD background 4j
      case 4: 
      eventfile="qcd_madgraph_4j_14tev_100k_gcuts.lhe";
      samplename="QCD4j";
      signal = false;
      xsec = .19922E+08; // pb
      xsec *= 1e3; // fb
      break;

    // ttbar productiom
    // Only fully hadronic decays considered
    // leptonic decays can be vetoed using leptons
      case 5:
      eventfile="qcd_madgraph_tt_hadr_14tev_100k_gcuts.lhe";
      samplename="QCDttbar";
      signal = false;
      xsec = .95450E+02; // pb
      xsec *= 1e3; // fb
      xsec *= 2.5; // NNLO K-factor
      break;


      default:
      std::cout<<"Invalid Monte Carlo sample, exit"<<std::endl;
      exit(-10);
      break;
    }

    // Initialize Pythia8
    eventfile= samples_path+eventfile;
    Pythia pythiaRun(std::string(PYTHIADIR)+"/xmldoc/");
    InitPythia(pythiaRun, eventfile);

    // Initialse Analyses for sample
    vector<Analysis*> sampleAnalyses;
    sampleAnalyses.push_back(new UCLAnalysis(samplename));
    sampleAnalyses.push_back(new DurhamAnalysis(samplename));
    
    int nev_tot = 0;

    // Begin loop over events
    for (int iEvent = 0; ;  ++iEvent) {

      nev_tot++;
      // Uncomment if prefer to run over a subset of events only
      //if(nev_tot>5e3) break;

      if (!pythiaRun.next()) {
      	// Stop showering when the end of the LHE file is reached
      	if (pythiaRun.info.atEndOfFile()) {
         cout << "Info: end of input file reached" << endl;
         break;
       }
     }

      // Obtain the final state
     finalState fs; get_final_state_particles(pythiaRun, fs);

     // Total analyses
     for (size_t i=0; i<HH4bAnalyses.size(); i++)
      HH4bAnalyses[i]->Analyse(signal, fs);

      // Sample analyses
     for (size_t i=0; i<sampleAnalyses.size(); i++)
      sampleAnalyses[i]->Analyse(signal, fs);

  }


  // Compute the total weight of the sample
    // Should coincide with nev_pass for btag prob of 1.0 and light jet mistag prob of 0.0
    for (size_t i=0; i<sampleAnalyses.size(); i++)
    {
      const double total_weight = sampleAnalyses[i]->GetWeight();
      const int nev_pass = sampleAnalyses[i]->GetNPassed();

      // Save results for cross-sections and number of events
      // Use LHC Run II and HL-LHC luminosities
      out_results<<"\nSample = "<< samplename<<" , Analysis = "<< sampleAnalyses[i]->GetName()<<std::endl;
      out_results<<"nev_tot(MC), nev_pass(MC) = "<<nev_tot<<" , "<<nev_pass<<std::endl;
      out_results<<"xsec_tot, xsec_pass (fb) = "<<xsec<< " , "<<total_weight/nev_tot<<std::endl;
      // LHC run II numbers
      out_results<<"nev_tot, nev_pass (300 1/fb) = "<< lumi_run2*xsec*total_weight/nev_tot<<std::endl;
      // HL-LHC numbers
      out_results<<"nev_tot, nev_pass (3000 1/fb) = "<< lumi_hllhc*xsec*total_weight/nev_tot<<std::endl;
  }

  // Free sample analyses
 for (size_t i=0; i<sampleAnalyses.size(); i++)
  delete sampleAnalyses[i];


}

// Finish up
out_results.close();
for (size_t i=0; i<HH4bAnalyses.size(); i++)
  delete HH4bAnalyses[i];
  
  // End of the main progream
return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
