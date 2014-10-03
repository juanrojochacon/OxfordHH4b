/////////////////////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>

//ROOT
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"
#include "TLatex.h"

// FastJet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
// FastJet contrib
#include "fastjet/contrib/VariableRPlugin.hh"

// Pythia8
#include "Pythia8/Pythia.h"

// namespaces
using namespace Pythia8; 
using namespace std;
using namespace fastjet;
using namespace contrib;

// General settings for analysis and jet reconstruction
#include "settings.h"

// Histogram plotting
#include "histograms.h"

// Analysis code
#include "analysis.h"

////////////////////////////////////////////////////////////////////////////////////////

int main() {
  
  cout<<"\n ***************************************************"<<endl;
  cout<<"\n Double Higgs Production in the 4b final state \n "<<endl;
  cout<<"***************************************************\n"<<endl;
  
  // Output some of the settings
  cout<<"\n ***************************************************\n"<<endl;
  cout<<"jetR (small R)= "<<jetR_0p5<<std::endl;
  cout<<"jetR (large R) = "<<jetR_1p2<<std::endl;
  cout<<"btag_prob = "<<btag_prob<<std::endl;
  cout<<"btag_mistag = "<<btag_mistag<<std::endl;
  cout<<"\n ***************************************************"<<endl;

  // Set here the path to the MC samples Dropbox folder
  string samples_path="/Users/juanrojo/Dropbox/HH4bMC/";

  // Select the analysis strategy
  // Here one can add more analysis strategies keeping the same I/O format
  cout<<"\n ***************************************************\n"<<endl;
  cout<<"Please select the analysis strategy"<<std::endl;
  cout<<"Enter 1 for UCL cut-based "<<std::endl;
  cout<<"Enter 2 for Durham substructure-based "<<std::endl;
  cout<<"\n ***************************************************"<<endl;
  int strat;
  cin>>strat;
  string strategy;
  if(strat==1)strategy="ucl";
  if(strat==2)strategy="durham";

  // Init histograms for the plots
  histo_init();

  // Open file to output the results
  ofstream out_results;
  out_results.open("hh4b.res");
  
  /* ---------------------------------------------------------------------------
  //
  // Loop over signal and background events
  // Signal: gg -> hh with mg5_amc at LO in EFT with heavy quark mass effects
  // Background: QCD 4b, 2b2j, ttbar fully hadronic
  // Need to add new sources of background, in particular QCD 4j
  // Other backgrounds, like single Higgs, also need to be included
  //
  ---------------------------------------------------------------------------*/
  
  int const nlhe=5; // Number of signal and bkg MC samples

  for(int ilhe=0; ilhe<nlhe;ilhe++){
    
    // Create histograms for this particular sample
    histo_create();
    
    // Name of the input file
    string eventfile;
    
    // Now we loop over the MC samples
    // The total cross-section for each sample can be
    // read from the .lhe file, see the number next to:
    // #  Integrated weight (pb)  
    double xsec = 0.9;

    // SM gg->HH, 100K
    if(ilhe==0) {
      eventfile="HH_sm_eft_100K.lhe";
      xsec = 0.17291e-1;// pb
      xsec *= 1e3; // fb
      xsec *= 0.34; // HH-> 4b BR
      xsec *= 2.5; // NNLO K-factor
    }

    // BSM gg->HH, 100K, trilinear multiplied by a factor 10
    else if(ilhe==1) {
      eventfile="HH_bsm_lam_10_eft_100K.lhe";
      xsec = .27923E+00;// pb
      xsec *= 1e3; // fb
      xsec *= 0.34; // HH-> 4b BR
      xsec *= 2.5; // NNLO K-factor
    }

    // QCD background, 4b
    else if(ilhe==2) {
      eventfile="qcd_madgraph_4b_14tev_100k_gcuts.lhe";
      xsec = .58013e+03; // pb
      xsec *= 1e3; // fb
      // Need to add NLO K-factor here
    }

    // QCD background, 2b2j
    else if(ilhe==3) {
      eventfile="qcd_madgraph_2b2j_14tev_100k_gcuts.lhe";
      xsec = .53822E+06; // pb
      xsec *= 1e3; // fb
      // Need to add NLO K-factor here
    }

    // ttbar productiom
    // Only fully hadronic decays considered
    // leptonic decays can be vetoed using leptons
    else if(ilhe==4) {
      eventfile="qcd_madgraph_tt_hadr_14tev_100k_gcuts.lhe";
      xsec = .95450E+02; // pb
      xsec *= 1e3; // fb
      xsec *= 2.5; // NNLO K-factor
    }
    
    else{
      std::cout<<"Invalid Monte Carlo sample, exit"<<std::endl;
      exit(-10);
    }

    // Write some info
    eventfile= samples_path+eventfile;
    cout<<"\n **********************************************************\n"<<endl;
    cout<<"Sample = "<<eventfile<<std::endl;
    cout<<"xsec (fb) = "<<xsec<<std::endl;
    cout<<"\n ***********************************************************"<<endl;

    // Initialize Pythia8
    Pythia pythiaRun;
    InitPythia(pythiaRun, eventfile);
    
    // Counters
    int nev_tot=0, nev_pass=0;
    // event weights - account for b tag and light jet mistag probabilities
    vector<double> weights;
    
    // Begin loop over events
    for (int iEvent = 0; ;  ++iEvent) {
      
      nev_tot++;
      // Uncomment if prefer to run over a subset of events only
      // if(nev_tot>5e3) break;

      if (!pythiaRun.next()) {
	// Stop showering when the end of the LHE file is reached
	if (pythiaRun.info.atEndOfFile()) {
	  cout << "Info: end of input file reached" << endl;
	  break;
	}
      }
      
      // Perform the jet clustering
      // Also b tagging and substructure done here
      // Returns the weight of the event including b tag and mistag probabilities
      // Each analysis strategy used a specific jet clustering step
      vector<fastjet::PseudoJet>  bjets;
      vector<PseudoJet> higgs_candidates;
      double event_weight=0; // initialization

      // UCL strategy - cut based, small R jets
      if(strategy=="ucl"){
	jet_clustering_analysis_smallR(pythiaRun, bjets, event_weight );
      }
      // Durham strategy - substructure, large R jets
      else if(strategy=="durham"){
	jet_clustering_analysis_largeR(pythiaRun, higgs_candidates, event_weight );
      }
      else{
	std::cout<<"Invalid analysis strategy"<<std::endl;
	exit(-10);
      }
      
      // Skip this event if basic selection fails, that is, if
      // the event weight vanishes
      if(event_weight<1e-30) continue;

      // Now here comes the second piece of analysed based on the reconstructed jets
      bool tagging = false;

      // UCL strategy - cut based, small R jets
      if(strategy=="ucl"){
	tagging = analysis_4b_ucl(bjets, event_weight);
      }
      // Durham strategy - substructure, large R jets
      else if(strategy=="durham"){
	tagging = analysis_4b_durham(higgs_candidates, event_weight);
      }
      else{
	std::cout<<"Invalid analysis strategy"<<std::endl;
	exit(-10);
      }

      if(tagging){
	// This event has passed all cuts
	// and has been tagged as HH->4b, so we keep its weight
	nev_pass++;
	weights.push_back(event_weight);
      }
      
    }

    // Check consistency in number of events that pass
    if(nev_pass != weights.size()) exit(-10);
    
    // Compute the total weight of the sample
    // Should coincide with nev_pass for btag prob of 1.0 and light jet mistag prob of 0.0
    double total_weight=0;
    for(unsigned i=0; i < weights.size(); i++) total_weight += weights.at(i);
    std::cout<<"nev_pass, total_weight = "<<nev_pass<<"\t "<<total_weight<<std::endl;

    // Save results for cross-sections and number of events
    // Use LHC Run II and HL-LHC luminosities
    out_results<<"\nSample = "<< eventfile<<std::endl;
    out_results<<"nev_tot(MC), nev_pass(MC) = "<<nev_tot<<" , "<<nev_pass<<std::endl;
    out_results<<"xsec_tot, xsec_pass (fb) = "<<xsec<< " , "<<total_weight/nev_tot<<std::endl;
    // LHC run II numbers
    out_results<<"nev_tot, nev_pass (300 1/fb) = "<< lumi_run2*xsec*total_weight/nev_tot<<std::endl;
    // HL-LHC numbers
    out_results<<"nev_tot, nev_pass (3000 1/fb) = "<< lumi_hllhc*xsec*total_weight/nev_tot<<std::endl;

    // Create the histograms for this sample
    histo_plot();
    
  } // End loop over signal and background samples
  
    // Plot the final combined histograms
  histo_plot_final(); 
  
  // End of the main progream
  return 0;
  
}
  
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
