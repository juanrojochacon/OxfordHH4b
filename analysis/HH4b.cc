//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

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

// Pythia8
#include "Pythia8/Pythia.h"

// namespaces
using namespace Pythia8; 
using namespace std;
using namespace fastjet;

// General settings for analysis
// and jet reconstruction
#include "settings.h"

// Routines for b and tau tagging
#include "tagging.h"

// Histograms
#include "histograms.h"

// Analysis
#include "analysis.h"


///////////////////////////////////////////////////////////

void InitPythia(Pythia & pythiaRun, string eventfile){

  // Initialize random seed
  srand (time(NULL));
  std::cout<<"time = "<<time(NULL)<<std::endl;
  double random = double(rand())/RAND_MAX;
  std::cout<<"\n\n Random number I = "<<random<<"\n\n"<<std::endl;
  random = double(rand())/RAND_MAX;
  std::cout<<"\n\n Random number II = "<<random<<"\n\n"<<std::endl;
  //  exit(-10);

  // Random seed
  pythiaRun.readString("Random:setSeed = on");
  double random_seed_pythia = 100000 * double(rand())/RAND_MAX;
  ostringstream o;
  o<<"Random:seed = "<<int(random_seed_pythia);
  cout<<o.str()<<endl;
  pythiaRun.readString(o.str());

  // Initialize Les Houches Event File run. List initialization information.
  pythiaRun.readString("Beams:frameType = 4"); 
  
  // Shower settings
  pythiaRun.readString("SpaceShower:QEDshowerByQ  = off"); // QED shower offf
  pythiaRun.readString("SpaceShower:QEDshowerByL  = off"); // QED shower offf
  pythiaRun.readString("HadronLevel:all = off"); // Of hadronization
  pythiaRun.readString("TimeShower:QEDshowerByQ = off");  // QED off on ISR / quarks irradiate photons
  pythiaRun.readString("TimeShower:QEDshowerByL = off");  // QED off on ISR / leptons irradiate photons  
  pythiaRun.readString("TimeShower:QEDshowerByGamma = off");  // Allow photons to branch into lepton or quark pairs 
  pythiaRun.readString("PartonLevel:MI = off"); // Off multiple interactions (UE) 
  pythiaRun.readString("PartonLevel:ISR = on");  // Shower on
  pythiaRun.readString("PartonLevel:FSR = on");  // Shower on

  // Higgs decays always into 4b
  pythiaRun.readString("25:onMode = off");
  pythiaRun.readString("25:onIfAll = 5 -5");

  // b quarks and do not decay
  pythiaRun.readString("5:mayDecay = no");
  pythiaRun.readString("-5:mayDecay = no");

  // Read the Les Houches Event File
  string ofile;
  ofile="Beams:LHEF = Events/"+eventfile;
  pythiaRun.readString(ofile.c_str());

   // Initialization
  pythiaRun.init();

}
/////////////////////////////////////////////////////////////

 
int main() {

  cout<<"\n Double Higgs Production in the 4b final state \n "<<endl;


  // Init histogarms
  histo_init();

  // Open file to output the results
  ofstream out_results;
  out_results.open("hh4b.res");

  // Loop over signal and background events
  // Signal: gg -> hh with mg5_amc at LO in EFT with heavy quark mass effects
  // Background: QCD 4b, 2b2j, 4j production
  // Both signal and backgrounds have been matched to parton showers
  int const nlhe=3;
  for(int ilhe=0; ilhe<nlhe;ilhe++){

    // Name of the input file
    string eventfile;
    double xsec;

    // Create histograms
    histo_create();

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
      eventfile="HH_BMS_Lam_10_eft_100K_gcut.lhe";
      xsec = .27923E+00;// pb
      xsec *= 1e3; // fb
      xsec *= 0.34; // HH-> 4b BR
      xsec *= 2.5; // NNLO K-factor
    }
    else if(ilhe==2) {
      eventfile="qcd_madgraph_4b_14tev_100k_gcuts.lhe";
      xsec = .58013e+03; // pb
      xsec *= 1e3; // fb
    }
    else if(ilhe==3) {
      eventfile="qcd_madgraph_2b2j_14tev_100k_gcuts.lhe";
      xsec = .53822E+06; // pb
      xsec *= 1e3; // fb
    }
    else{
      std::cout<<"Invalid value of LHEF"<<std::endl;
      exit(-10);
    }
    
    // Init Pythia8
    Pythia pythiaRun;
    InitPythia(pythiaRun, eventfile);

    // Counter
    int nev_pass=0;
    int nev_pass_boosted=0;
    int nev_tot=0;
    
    // Begin loop over events
    for (int iEvent = 0; ;  ++iEvent) {
      
      nev_tot++;

      if(iEvent > 1e5) break;

      if (!pythiaRun.next()) {
	if (pythiaRun.info.atEndOfFile()) {
	  cout << "Info: end of input file reached" << endl;
	  if(ilhe<3)break;
	  if(ilhe==3){
	    // Begin again from the same file
	    // but different random seed
	    Pythia pythiaRun;
	    std::cout<<"Using again file with different random seed = "<<eventfile<<std::endl;
	    InitPythia(pythiaRun, eventfile);
	  }
	}
      }
      
      vector<PseudoJet> particles;
      
      // Read event
      double px = 0;
      double py = 0;
      double pz =0;
      double E = 0;
      
      for (int i = 0; i < pythiaRun.event.size(); i++){
	int particle_id = pythiaRun.event[i].id();
	int particle_status = pythiaRun.event[i].status();
	
	//      std::cout<<particle_id<<" "<<particle_status<<std::endl;
	
	// Check final state particles
	if( particle_status > 0 ) {
	  
	  // Get the particle kinematics
	  px= pythiaRun.event[i].px();
	  py= pythiaRun.event[i].py();
	  pz= pythiaRun.event[i].pz();
	  E= pythiaRun.event[i].e();
	  
	  // quarks and gluons
	  // including b-quarks
	  if(abs(particle_id)<6 || particle_id==21 ){
	    particles.push_back( fastjet::PseudoJet(px,py,pz,E) );
	    particles.at(particles.size()-1).set_user_index(particle_id);
	  }
	  // beam remnants
	  else if(particle_id > 2000 ){
	    particles.push_back( fastjet::PseudoJet(px,py,pz,E) );
	    particles.at(particles.size()-1).set_user_index(particle_id);
	  }
	  else{
	    std::cout<<"Invalid particle ID = "<<particle_id<<std::endl;
	    exit(-10);
	  }
	  
	} // End loop over final state particles
	
      } // End loop over particles in event
      
      // Check event energy conservation here
      double px_tot=0;
      double py_tot=0;
      double pz_tot=0;
      double E_tot=0;
      
      // quarks and gluons
      for(unsigned ij=0;ij<particles.size();ij++){
	px_tot+= particles.at(ij).px();
	py_tot+= particles.at(ij).py();
	pz_tot+= particles.at(ij).pz();
	E_tot+= particles.at(ij).E();
      }
      // Check energy-momentum conservation
      double const Eref=14000;
      double const tol_emom=1.0;
      if( fabs(px_tot) > tol_emom || fabs(py_tot)  > tol_emom || fabs(pz_tot)  > tol_emom || fabs(E_tot-Eref)  > tol_emom ){
	std::cout<<"\n ********************************************************************** \n"<<std::endl;
	std::cout<<"No conservation of energy in Pythia after shower "<<std::endl;
	std::cout<<"px_tot = "<<px_tot<<std::endl;
	std::cout<<"py_tot = "<<py_tot<<std::endl;
	std::cout<<"pz_tot = "<<pz_tot<<std::endl;
	std::cout<<"E_tot, Eref = "<<E_tot<<" "<<Eref<<std::endl;
	exit(-10);
	std::cout<<"\n ********************************************************************** \n"<<std::endl;
      }
      
      // Do some simple jet clustering with FastJet
      // Jet clustering
      JetDefinition akt(antikt_algorithm, jetR);
      // Cluster all particles
      ClusterSequence cs_akt(particles, akt);
      // Get all the jets (no pt cut here)
      vector<fastjet::PseudoJet> jets_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
      
      // Check again four-momentum conservation, this time applied to jets
      // quarks and gluons
      px_tot=0;
      py_tot=0;
      pz_tot=0;
      E_tot=0;
      for(unsigned ij=0;ij<jets_akt.size();ij++){
	px_tot+= jets_akt.at(ij).px();
	py_tot+= jets_akt.at(ij).py();
	pz_tot+= jets_akt.at(ij).pz();
	E_tot+= jets_akt.at(ij).E();
      }
      
      // Check energy-momentum conservation
      if( fabs(px_tot) > tol_emom || fabs(py_tot)  > tol_emom || fabs(pz_tot)  > tol_emom || fabs(E_tot-Eref)  > tol_emom ){
	std::cout<<"\n ********************************************************************** \n"<<std::endl;
	std::cout<<"No conservation of energy in Pythia after shower and jet reconstruction "<<std::endl;
	std::cout<<"px_tot = "<<px_tot<<std::endl;
	std::cout<<"py_tot = "<<py_tot<<std::endl;
	std::cout<<"pz_tot = "<<pz_tot<<std::endl;
	std::cout<<"E_tot, Eref = "<<E_tot<<" "<<Eref<<std::endl;
	exit(-10);
	std::cout<<"\n ********************************************************************** \n"<<std::endl;
      }

      // Insert here your analysis
      // based on the reconstructed jets 
      // b-tagging effects are included here
      // Returns true if the event is tagged as arising from HH->4b
      bool hh4b_tagged = bbbb_analysis(jets_akt);
      if(hh4b_tagged) nev_pass++;

      // Here second analysis
      // Based on jet substructure
      // Avoid double counting as compared to the resolved regime
      if(!hh4b_tagged){
	bool hh4b_tagged_boosted = bbbb_analysis_boosted(jets_akt);
	if(hh4b_tagged_boosted) nev_pass_boosted++;
      }      

      // Output intermediate results every tot generations
      if( ( nev_tot % 10000 ) == 0 ){
	std::cout<<"\n\nSample = "<< eventfile<<std::endl;
	std::cout<<"Resolved: "<<std::endl;
	std::cout<<"nev_tot(MC), nev_pass(MC) = "<<nev_tot<<" , "<<nev_pass<<std::endl;
	std::cout<<"Boosted: "<<std::endl;
	std::cout<<"nev_tot(MC), nev_pass(MC) = "<<nev_tot<<" , "<<nev_pass_boosted<<std::endl;
      }         

    } // End loop over events in a given sample
 
    double const lumi=300;
   
    out_results<<"\n\nSample = "<< eventfile<<std::endl;
    out_results<<"\n Resolved: "<<std::endl;
    out_results<<"nev_tot(MC), nev_pass(MC) = "<<nev_tot<<" , "<<nev_pass<<std::endl;
    out_results<<"xsec_tot, xsec_pass (fb) = "<<xsec<< " , "<<xsec*double(nev_pass)/nev_tot<<std::endl;
    out_results<<"nev_tot(300 1/fb), nev_pass (300 1/fb) = "<<xsec*lumi<< " , "<<lumi*xsec*double(nev_pass)/nev_tot<<std::endl;
    out_results<<"nev_tot(3000 1/fb), nev_pass (3000 1/fb) = "<<xsec*lumi*10<< " , "<<lumi*10*xsec*double(nev_pass)/nev_tot<<std::endl;

    out_results<<"\n Boosted: "<<std::endl;
    out_results<<"nev_tot(MC), nev_pass(MC) = "<<nev_tot<<" , "<<nev_pass_boosted<<std::endl;
    out_results<<"xsec_tot, xsec_pass (fb) = "<<xsec<< " , "<<xsec*double(nev_pass_boosted)/nev_tot<<std::endl;
    out_results<<"nev_tot(300 1/fb), nev_pass (300 1/fb) = "<<xsec*lumi<< " , "<<lumi*xsec*double(nev_pass_boosted)/nev_tot<<std::endl;
    out_results<<"nev_tot(3000 1/fb), nev_pass (3000 1/fb) = "<<xsec*lumi*10<< " , "<<lumi*10*xsec*double(nev_pass_boosted)/nev_tot<<std::endl;

    double scalefactor=1.0;
    histo_plot(scalefactor);

  } // End loop over signal and background samples
  
  // Very final combined joint plot
  histo_plot_final(); 
  
  // End of the main progream
  return 0;
  
}
  
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
