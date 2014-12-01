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

#include "ucl.h"
#include "ucl_vr.h"
#include "durham.h"
#include "oxford_res_vr.h"
#include "oxford_res_fr.h"
#include "oxford_boost_vr.h"
#include "oxford_boost_fr.h"
#include "variableR.h"
#include "amcatnlo.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"


using namespace Pythia8;

////////////////////////////////////////////////////////////////////////////////////////

int main() 
{  
// Results output
  ofstream out_results;
  const string outDir = "./" + std::string(RESDIR) + "/hh4b.res";
  out_results.open(outDir.c_str());

// Set here the path to the MC samples Dropbox folder
  string samples_path=std::string(SAMPLEDIR);

// Initialise analyses
  vector<Analysis*> HH4bAnalyses;
  vector<Analysis*> signalAnalyses;
  vector<Analysis*> backgroundAnalyses;

// aMC@NLO numbers for total xsecs
  HH4bAnalyses.push_back(new AMCAnalysis("total"));

  HH4bAnalyses.push_back(new UCLAnalysis("total"));
  HH4bAnalyses.push_back(new DurhamAnalysis("total"));
  HH4bAnalyses.push_back(new UCLVRAnalysis("total"));
  HH4bAnalyses.push_back(new OxfordResVRAnalysis("total"));
  HH4bAnalyses.push_back(new OxfordResFRAnalysis("total"));
  HH4bAnalyses.push_back(new OxfordBoostVRAnalysis("total"));
  HH4bAnalyses.push_back(new OxfordBoostFRAnalysis("total"));

  signalAnalyses.push_back(new UCLAnalysis("signal"));
  signalAnalyses.push_back(new DurhamAnalysis("signal"));
  signalAnalyses.push_back(new UCLVRAnalysis("signal"));
  signalAnalyses.push_back(new OxfordResVRAnalysis("signal"));
  signalAnalyses.push_back(new OxfordResFRAnalysis("signal"));
  signalAnalyses.push_back(new OxfordBoostVRAnalysis("signal"));
  signalAnalyses.push_back(new OxfordBoostFRAnalysis("signal"));

  backgroundAnalyses.push_back(new UCLAnalysis("background"));
  backgroundAnalyses.push_back(new DurhamAnalysis("background"));
  backgroundAnalyses.push_back(new UCLVRAnalysis("background"));
  backgroundAnalyses.push_back(new OxfordResVRAnalysis("background"));
  backgroundAnalyses.push_back(new OxfordResFRAnalysis("background"));
  backgroundAnalyses.push_back(new OxfordBoostVRAnalysis("background"));
  backgroundAnalyses.push_back(new OxfordBoostFRAnalysis("background"));



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

for(int ilhe=3; ilhe<nlhe;ilhe++){

// Name of the input file
  string eventfile;
  string samplename;
  bool signal = false;
  bool hepmc = false;

// Now we loop over the MC samples
// The total cross-section for each sample can be
// read from the .lhe file, see the number next to:
// #  Integrated weight (pb)  
double xsec_norm = 1e3; // pb -> fb conversion

// Number of events in sample
int nevt_sample = 0;

switch (ilhe)
{
// SM gg->HH, 100K
case 2: 
eventfile="HH_sm_eft_100K.lhe";
samplename="SMggHH";
signal = true;
hepmc = false;
nevt_sample = 1E5;
xsec_norm *= 0.34; // HH-> 4b BR
xsec_norm *= 2.5; // NNLO K-factor
break;

// BSM gg->HH, 100K, trilinear multiplied by a factor 10
case 3: 
eventfile="HH_bsm_lam_10_eft_100K.lhe";
samplename="BSMggHH";
signal = true;
hepmc = false;
nevt_sample = 1E5;
xsec_norm *= 0.34; // HH-> 4b BR
xsec_norm *= 2.5; // NNLO K-factor
continue;
break;

// QCD background, 4b
case 2: 
eventfile="mg5_qcd4b_14tev_combined.lhe";
samplename="QCD4b";
signal = false;
hepmc = false;
nevt_sample = 829431;
// Need to add NLO K-factor here
break;

// QCD background, 2b2j
case 3: 
eventfile="qcd_madgraph_2b2j_14tev_100k_gcuts.lhe";
samplename="QCD2b2j";
signal = false;
hepmc = false;
nevt_sample = 1060000;
// Need to add NLO K-factor here
break;

// QCD background 4j
case 5: 
eventfile="qcd_madgraph_4j_14tev_100k_gcuts.lhe";
samplename="QCD4j";
signal = false;
hepmc = false;
nevt_sample = 1E5;
break;

// ttbar productiom
// Only fully hadronic decays considered
// leptonic decays can be vetoed using leptons
case 6:
eventfile="qcd_madgraph_tt_hadr_14tev_100k_gcuts.lhe";
samplename="QCDttbar";
signal = false;
hepmc = false;
nevt_sample = 1E5;
xsec_norm *= 2.5; // NNLO K-factor
break;

// SM gg->HH, 100K
case 7: 
eventfile="2b2j.hepmc";
samplename="Sherpa2b2j";
signal = false;
hepmc = true;
nevt_sample = 5E4;
break;

// SM gg->HH, 100K
case 8: 
eventfile="4b.hepmc";
samplename="Sherpa4b";
signal = false;
hepmc = true;
nevt_sample = 4E4;
break;

default:
std::cout<<"Invalid Monte Carlo sample, exit"<<std::endl;
exit(-10);
break;
}

// Maximum number of events to loop over
const int nevt_max = min(nevt_sample,(int)5E4);

// Initialize inputs
eventfile= samples_path+eventfile;
std::cout << "Reading samples from: "<<eventfile<<std::endl;

Pythia pythiaRun(std::string(PYTHIADIR)); // Pythia input

// Initialise the event sample and unit weight
double unit_weight = 0;
if (!hepmc)
{
  InitPythia(pythiaRun, eventfile);

// Verify pythia event sample
  if (pythiaRun.info.nProcessesLHEF() != 1)
  {
    std::cerr << "Error: number of subprocesses in LHE sample is not equal to 1 " <<std::endl;
    std::cerr << "This code does not support this at the moment"<<std::endl;
    exit(-1);
  }
  
  // Event weight information from Pythia
  const double pythia_wgt = pythiaRun.info.sigmaLHEF(0); // Total sample weight
  const double xsec = xsec_norm*pythia_wgt;              // Normalised sample weight (fb, Kfactors)
  unit_weight = xsec/((double) nevt_max);      // Unit weight
} 
else
{
  std::ifstream hepmc_is( eventfile );            // HepMC input
  for (int iEvent = 0; iEvent < nevt_max; iEvent++)
      {
        if (! hepmc_is )
        {
          cerr << "Error: HepMC end of file!"<<endl;
          exit(-1);
        }

        HepMC::GenEvent evt;
        evt.read( hepmc_is );
        if( !evt.is_valid() ) 
        {
          cerr << "Error: Invalid HepMC event!" << endl;
          exit(-1);
        }

        // Update to generated xsec
        unit_weight = evt.cross_section()->cross_section();
      } 

      cout << "Event Generator xsec: "<< unit_weight*xsec_norm<<" (fb)"<< endl;
      unit_weight *= xsec_norm/nevt_max;

      // Reset istream
      hepmc_is.close();
}    

cout << "Unit Weight: "<< unit_weight<<" N_evt: "<<nevt_max<<endl;

// Initialse Analyses for sample
vector<Analysis*> sampleAnalyses;
sampleAnalyses.push_back(new AMCAnalysis(samplename));

sampleAnalyses.push_back(new UCLAnalysis(samplename));
sampleAnalyses.push_back(new DurhamAnalysis(samplename));
sampleAnalyses.push_back(new UCLVRAnalysis(samplename));
sampleAnalyses.push_back(new OxfordResVRAnalysis(samplename));
sampleAnalyses.push_back(new OxfordResFRAnalysis(samplename));
sampleAnalyses.push_back(new OxfordBoostVRAnalysis(samplename));
sampleAnalyses.push_back(new OxfordBoostFRAnalysis(samplename));

// counters
int nev_tot = 0;
double sample_xsec = 0;

double hepmc_weights = 0;
double hepmc_tries = 0;

// Begin loop over events
std::ifstream hepmc_is( eventfile );            // HepMC input
for (int iEvent = 0; ; ++iEvent) 
{
// Increment total counter and break if reached maximum
  nev_tot++; if(nev_tot>nevt_max) break;
  if (nev_tot % 1000 == 0 && hepmc) cout << "Event "<< nev_tot<< " processed." <<endl;
// Event information
  finalState fs;
  if (!hepmc)
  {
    if (!pythiaRun.next()) 
      if (pythiaRun.info.atEndOfFile()) // Stop showering when the end of the LHE file is reached
      {
        cout << "Info: end of input file reached" << endl;
        break;
      }

    // Obtain the final state
    get_final_state_particles(pythiaRun, fs);

  } 
  else  // HepMC format
  {
    HepMC::GenEvent evt;
    if ( hepmc_is ) 
    {
      evt.read( hepmc_is );
      if( !evt.is_valid() ) 
      {
        cerr << "Error: Invalid HepMC event!" << endl;
        exit(-1);
      }
    } else {
        cout << "Info: end of input file reached" << endl;
        break;
    }

    // Obtain the final state
    get_final_state_particles(evt, fs);
  }

// Calculate the sample cross-section
sample_xsec += unit_weight;

// Total analyses
for (size_t i=0; i<HH4bAnalyses.size(); i++)
  HH4bAnalyses[i]->Analyse(signal, unit_weight, fs);

// Sample analyses
for (size_t i=0; i<sampleAnalyses.size(); i++)
  sampleAnalyses[i]->Analyse(signal, unit_weight, fs);

// Signal
if (signal)
  for (size_t i=0; i<signalAnalyses.size(); i++)
    signalAnalyses[i]->Analyse(signal, unit_weight, fs);
  else
    for (size_t i=0; i<backgroundAnalyses.size(); i++)
      backgroundAnalyses[i]->Analyse(signal, unit_weight, fs);

  }

// Compute the total weight of the sample
// Should coincide with nev_pass for btag prob of 1.0 and light jet mistag prob of 0.0
  for (size_t i=0; i<sampleAnalyses.size(); i++)
  {
    const double wgt_pass = sampleAnalyses[i]->GetWeight();
    const int nev_pass = sampleAnalyses[i]->GetNPassed();

// Save results for cross-sections and number of events
// Use LHC Run II and HL-LHC luminosities
    out_results<<"\nSample = "<< samplename<<" , Analysis = "<< sampleAnalyses[i]->GetName()<<std::endl;
    out_results<<"nev_tot(MC), nev_pass(MC) = "<<nev_tot<<" , "<<nev_pass<<std::endl;
    out_results<<"xsec_tot, xsec_pass (fb) = "<<sample_xsec<< " , "<<wgt_pass<<std::endl;
// LHC run II numbers
    out_results<<"pass weight (300 1/fb) = "<< lumi_run2*wgt_pass<<std::endl;
// HL-LHC numbers
    out_results<<"pass weight, nev_pass (3000 1/fb) = "<< lumi_hllhc*wgt_pass<<std::endl;

// Save results for cross-sections and number of events
// Use LHC Run II and HL-LHC luminosities
    const double notCounted = sample_xsec - (sampleAnalyses[i]->GetCutWeight() + sampleAnalyses[i]->GetWeight());
    cout <<"\nSample = "<< samplename<<" , Analysis = "<< sampleAnalyses[i]->GetName()<<std::endl;
    cout <<"nev_tot(MC), nev_pass(MC) = "<<nev_tot<<" , "<<nev_pass<<std::endl;
    cout <<"xsec_tot, xsec_pass (fb) = "<<sample_xsec<< " , "<<wgt_pass<<std::endl;
    cout << "cutWeight: "<< sampleAnalyses[i]->GetCutWeight() <<" PassedWeight: "<<sampleAnalyses[i]->GetWeight()<<endl;
    cout << "Unaccounted for weights: "<< notCounted <<endl;

    if (notCounted > 1E-10)
    {
      cerr << "Warning: Unaccounted for weights are too large! "<<endl;
//exit(-1);
    }

    // Close stream
    hepmc_is.close();

  }


// Free sample analyses
  for (size_t i=0; i<sampleAnalyses.size(); i++)
    delete sampleAnalyses[i];

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
