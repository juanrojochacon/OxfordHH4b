// oxford_res_fr.cc

#include <algorithm>

#include "oxford_combined_rw.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
// #include "fastjet/contrib/VariableRPlugin.hh"

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

// pT cut for constituent b-quarks (GeV)
const double pt_btagging = 15;

OxfordCombinedRWAnalysis::OxfordCombinedRWAnalysis(std::string const& sampleName):
Analysis("oxford_combined_rw", sampleName)
{
  // ********************* Histogram settings******************

  const double DeltaRmin = 0;
  const double DeltaRmax = 5;

  const double DeltaPhimin = -3.2;
  const double DeltaPhimax = 3.2;

  const double DeltaEtamin = -2.5;
  const double DeltaEtamax = 2.5;
  
  const double m_min = 0.;
  const double m_max = 180.; 
 
  const double pt_min = 0.;
  const double pt_max = 900.;
  
  const double m_HH_min = 0.;
  const double m_HH_max = 600.; 
  
  const double pt_HH_min = 0.;
  const double pt_HH_max = 300.;
  
  const int nbins = 30;
  
  // ********************* Histogram definitions ******************

  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_res_C0");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_inter_C0");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_boost_C0");

  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_res_C0");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_inter_C0");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_boost_C0"); 
  
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_res_C0");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_inter_C0");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_boost_C0");

  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_res_C0");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_inter_C0");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_boost_C0"); 
  
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_res_C0");
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_inter_C0");
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_boost_C0"); 

  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_res_C0");
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_inter_C0");
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_boost_C0"); 

  // C1: Higgs histograms
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_res_C1");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_inter_C1");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_boost_C1");

  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_res_C1");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_inter_C1");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_boost_C1"); 
  
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_res_C1");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_inter_C1");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_boost_C1");

  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_res_C1");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_inter_C1");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_boost_C1"); 
  
  // C1: Di-Higgs histograms
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_res_C1");
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_inter_C1");
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_boost_C1"); 
  
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_res_C1");
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_inter_C1");
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_boost_C1"); 
  
  BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH_res_C1");
  BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH_inter_C1");
  BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH_boost_C1");

  BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH_res_C1");
  BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH_inter_C1");
  BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH_boost_C1");
  
  BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH_res_C1");
  BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH_inter_C1");
  BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH_boost_C1");

  
  // C2: Higgs histograms
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_res_C2");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_inter_C2");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_boost_C2");

  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_res_C2");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_inter_C2");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_boost_C2");
  
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_res_C2");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_inter_C2");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_boost_C2");

  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_res_C2");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_inter_C2");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_boost_C2"); 

  // C2: Di-Higgs histograms
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_res_C2");
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_inter_C2");
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_boost_C2");
  
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_res_C2");
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_inter_C2");
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_boost_C2"); 
  
  BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH_res_C2");
  BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH_inter_C2");
  BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH_boost_C2");

  BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH_res_C2");
  BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH_inter_C2");
  BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH_boost_C2");
  
  BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH_res_C2");
  BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH_inter_C2");
  BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH_boost_C2");
  
  
  // C3: Higgs histograms
  BookHistogram(new YODA::Histo2D(nbins, pt_min, pt_max, nbins, pt_min, pt_max), "ptHptH_res_C2");
  BookHistogram(new YODA::Histo2D(nbins, pt_min, pt_max, nbins, pt_min, pt_max), "ptHptH_inter_C2");
  BookHistogram(new YODA::Histo2D(nbins, pt_min, pt_max, nbins, pt_min, pt_max), "ptHptH_boost_C2");
  
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mHmH_res_C2");
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mHmH_inter_C2");
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mHmH_boost_C2");

  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_res_C3");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_inter_C3");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0_boost_C3");

  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_res_C3");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_inter_C3");
  BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1_boost_C3"); 
  
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_res_C3");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_inter_C3");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0_boost_C3");

  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_res_C3");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_inter_C3");
  BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1_boost_C3"); 
  
  // C3: Di-Higgs histograms
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_res_C3");
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_inter_C3");
  BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH_boost_C3");
  
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_res_C3");
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_inter_C3");
  BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH_boost_C3"); 
  
  BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH_res_C3");
  BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH_inter_C3");
  BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH_boost_C3");

  BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH_res_C3");
  BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH_inter_C3");
  BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH_boost_C3");
  
  BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH_res_C3");
  BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH_inter_C3");
  BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH_boost_C3");
  
  // 2D histograms
  BookHistogram(new YODA::Histo2D(nbins, pt_min, pt_max, nbins, pt_min, pt_max), "ptHptH_res_C3");
  BookHistogram(new YODA::Histo2D(nbins, pt_min, pt_max, nbins, pt_min, pt_max), "ptHptH_inter_C3");
  BookHistogram(new YODA::Histo2D(nbins, pt_min, pt_max, nbins, pt_min, pt_max), "ptHptH_boost_C3");
  
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mHmH_res_C3");
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mHmH_inter_C3");
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mHmH_boost_C3");
  
  // ********************* Cutflow histograms  **********************
  
  // Cross-section
  BookHistogram(new YODA::Histo1D( 4, 0, 4 ), "CF_res");
  BookHistogram(new YODA::Histo1D( 4, 0, 4 ), "CF_inter");
  BookHistogram(new YODA::Histo1D( 4, 0, 4 ), "CF_boost");

  // Number of events
  BookHistogram(new YODA::Histo1D( 4, 0, 4 ), "CFN_res");
  BookHistogram(new YODA::Histo1D( 4, 0, 4 ), "CFN_inter");
  BookHistogram(new YODA::Histo1D( 4, 0, 4 ), "CFN_boost");
  
  // ********************* Ntuple definition **********************
  const std::string tupleSpec = "# signal source weight";
  outputNTuple<<tupleSpec<<std::endl;
}

void OxfordCombinedRWAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
  // Only for signal for now
  //if (!signal) return;
  //if (fs.size() == 4) std::cout << "INFO: Apparently these are unshowered samples."<<std::endl; 

  Analysis::Analyse(signal, weightnorm, fs);

  // Set initial weight
  const double event_weight = weightnorm;

  
  // ********************************* Simple event categorisation  ***************************************
  //===============
  // Small-R jets
  //===============
  static double const ResJetR=0.4;
  fastjet::JetDefinition akt_res(fastjet::antikt_algorithm, ResJetR);
  fastjet::ClusterSequence cs_akt_res(fs, akt_res);
  std::vector<fastjet::PseudoJet> smallRJets = sorted_by_pt( cs_akt_res.inclusive_jets()  ); // Get all the jets (no pt cut here)
  
  // Basic kinematic cuts
  std::vector<fastjet::PseudoJet> smallRJetsSel;
  for( size_t i = 0; i < smallRJets.size(); i++){
    
    if( smallRJets.at(i).pt() < 40. ) continue;
    if( fabs( smallRJets.at(i).eta() ) > 2.5 ) continue;
    
    smallRJetsSel.push_back( smallRJets.at(i) );
  }
  

  smallRJetsSel = sorted_by_pt( smallRJetsSel  ); // Resort
  std::vector<bool> isFakeSR_vec;         // Vector specifying if each jet is fake or not
  BTagging( smallRJetsSel, isFakeSR_vec );

  if (smallRJetsSel.size() != isFakeSR_vec.size())
    std::cerr << "Error: vector sizes do not match! "<<std::endl;

  // Form small-R b-jets vector
  std::vector<fastjet::PseudoJet> bJets;
  const int nBTaggedJets = std::min(smallRJetsSel.size(), (size_t)4);   // Both for now
  const int nJets = std::min(smallRJetsSel.size(), (size_t)4);          // Both for now
  for( int i = 0; i < nBTaggedJets; i++)
    bJets.push_back(smallRJetsSel[i]);
  
  //=======================================================
  // Small-R track jets from charged final state particles
  //=======================================================
  finalState fsc;
  //Select only charged fs particles
  for(int i=0; i<(int)fs.size(); i++){
      int userid = fs.at(i).user_index();
      if(abs(userid)<6) fsc.push_back(fs.at(i));
  }
  
  double jetR = 0.3;
  fastjet::JetDefinition jd_subjets(fastjet::antikt_algorithm, jetR);
  fastjet::ClusterSequence cs_subjets(fsc, jd_subjets);
  std::vector<fastjet::PseudoJet> trackjets = sorted_by_pt( cs_subjets.inclusive_jets()  );
  
  // Basic kinematic cuts
  std::vector<fastjet::PseudoJet> trackjetsSel;
  for( size_t i = 0; i < trackjets.size(); i++){
    
    if( trackjets.at(i).pt() < 12 ) continue;
    if( fabs( trackjets.at(i).eta() ) > 2.5 ) continue;
    
    trackjetsSel.push_back( trackjets.at(i) );
  }
      
  //===============
  // Large-R jets
  //===============
  static double const BoostJetR=1.0;
  fastjet::JetDefinition akt_boost(fastjet::antikt_algorithm, BoostJetR);
  fastjet::ClusterSequence cs_akt_boost(fs, akt_boost);
  std::vector<fastjet::PseudoJet> largeRJets = sorted_by_pt( cs_akt_boost.inclusive_jets()  ); // Get all the jets (no pt cut here)
  std::vector<fastjet::PseudoJet> bbFatJets;
  
  // Basic kinematic cuts
  std::vector<fastjet::PseudoJet> largeRJetsSel;
  for( size_t i = 0; i < largeRJets.size(); i++){
 
    //-----------------------------------------------------
    // pt acceptance cut
    if( largeRJets.at(i).pt() < 200. ) continue;
    //-----------------------------------------------------
    // Pseudorapidity acceptance cut
    if( fabs( largeRJets.at(i).eta() ) > 2.0 ) continue;

    //-----------------------------------------------------
    // Check if jet is mass-drop tagged
    fastjet::JetDefinition CA10(fastjet::cambridge_algorithm, 1.0);
    fastjet::ClusterSequence cs_sub( largeRJets.at(i).constituents(), CA10);
    fastjet::PseudoJet ca_jet = sorted_by_pt(cs_sub.inclusive_jets())[0];

    // parameters are specified in settings.h
    fastjet::MassDropTagger md_tagger(mu, ycut);
    fastjet::PseudoJet tagged_jet = md_tagger(ca_jet);
    
    // Skip jet if it fails mass-drop tagging
    if (tagged_jet == 0 ) continue;

    //-----------------------------------------------------
 
    largeRJetsSel.push_back( largeRJets.at(i) );
  }

  largeRJetsSel = sorted_by_pt(largeRJetsSel);  // Resort
  std::vector<int> nBSubjetsLR_vec; // Vector specifying how many real b subjets there are
  BTaggingFJ( largeRJetsSel, trackjetsSel, nBSubjetsLR_vec);

  if( largeRJetsSel.size() != nBSubjetsLR_vec.size() )
    std::cout << "ERROR: b-tagging vector sizes don't match number of fat jets" << std::endl;
  


  const int nFatJets = (int)largeRJetsSel.size();
  const int nBBTaggedFatJets = std::min(largeRJetsSel.size(), (size_t)2); ;

  for( int i = 0; i < nBBTaggedFatJets; i++)
    bbFatJets.push_back( largeRJetsSel.at(i) );
  
  
  //=============================
  // C0: Cutflow before cuts
  // ============================
  FillHistogram("CF_res", event_weight, 0.1);
  FillHistogram("CF_inter", event_weight, 0.1);
  FillHistogram("CF_boost", event_weight, 0.1);

  FillHistogram("CFN_res", 1., 0.1);
  FillHistogram("CFN_inter", 1., 0.1);
  FillHistogram("CFN_boost", 1., 0.1);
  
  //===================================================
  // C1: Basic kinematic cuts + Higgs mass window cut
  //===================================================
  if( nFatJets >= 2 ){      
    
      // Higgs mass window cut
      if( fabs(largeRJetsSel[0].m() - 125.) < 40.0 && fabs(largeRJetsSel[1].m() - 125.) < 40.0 ){
	
	    // Record cutflow
	    FillHistogram("CF_boost", event_weight, 1.1);
	    FillHistogram("CFN_boost", 1., 1.1);
	    
	    // Histograms for reconstructed Higgs candidates
	    FillHistogram("pt_H0_boost_C1", event_weight, largeRJetsSel[0].pt());
	    FillHistogram("pt_H1_boost_C1", event_weight, largeRJetsSel[1].pt());
	    
	    FillHistogram("m_H0_boost_C1", event_weight, largeRJetsSel[0].m());
	    FillHistogram("m_H1_boost_C1", event_weight, largeRJetsSel[1].m());
	    
	    // Reconstruct di-Higgs system
	    const fastjet::PseudoJet dihiggs_boost = largeRJetsSel[0] + largeRJetsSel[1];
	    
	    // Histograms for reconstructed di-Higgs system
	    FillHistogram("m_HH_boost_C1", event_weight, dihiggs_boost.m());
	    FillHistogram("pt_HH_boost_C1", event_weight, dihiggs_boost.pt());
	    FillHistogram("dR_HH_boost_C1", event_weight, largeRJetsSel[0].delta_R(largeRJetsSel[1]) );
	    FillHistogram("dPhi_HH_boost_C1", event_weight, getDPhi(largeRJetsSel[0].phi(), largeRJetsSel[1].phi()) );
	    FillHistogram("dEta_HH_boost_C1", event_weight, fabs( largeRJetsSel[0].eta() - largeRJetsSel[1].eta()) );
	    

      }
  }
  if( nJets >= 4 ){
    
      // Reconstruct Higgs candidates from small-R jets
      std::vector<fastjet::PseudoJet> higgs_res;
      Reco_Resolved( smallRJetsSel, higgs_res );
    
      // Higgs mass window cut
      if( fabs(higgs_res[0].m() - 125.) < 40.0 && fabs(higgs_res[1].m() - 125.) < 40.0 ){
	
	  // Record cutflow
	  FillHistogram("CF_res", event_weight, 1.1);
	  FillHistogram("CFN_res", 1., 1.1);
	  
	  // Histograms for reconstructed Higgs candidates
	  FillHistogram("pt_H0_res_C1", event_weight, higgs_res[0].pt());
	  FillHistogram("pt_H1_res_C1", event_weight, higgs_res[1].pt());
	  
	  FillHistogram("m_H0_res_C1", event_weight, higgs_res[0].m());
	  FillHistogram("m_H1_res_C1", event_weight, higgs_res[1].m());
      
	  // Reconstruct di-Higgs system
	  const fastjet::PseudoJet dihiggs_res = higgs_res[0] + higgs_res[1];
	  
	  // Histograms for reconstructed di-Higgs system
	  FillHistogram("m_HH_res_C1", event_weight, dihiggs_res.m());
	  FillHistogram("pt_HH_res_C1", event_weight, dihiggs_res.pt());
	  FillHistogram("dR_HH_res_C1", event_weight, higgs_res[0].delta_R(higgs_res[1]) );
	  FillHistogram("dPhi_HH_res_C1", event_weight, getDPhi(higgs_res[0].phi(), higgs_res[1].phi()) );
	  FillHistogram("dEta_HH_res_C1", event_weight, fabs( higgs_res[0].eta() - higgs_res[1].eta()) );
      }
  }
  if( nJets >= 2 &&  nFatJets == 1 ){
    
    // Reconstruct Higgs candidates from large-R and small-R jets
    std::vector<fastjet::PseudoJet> higgs_inter;
    int nBJets_SR = 0;
    bool isRecoInter = Reco_Intermediate( smallRJetsSel, isFakeSR_vec, largeRJetsSel[0], nBJets_SR, higgs_inter );

    // Check if reconstruction was successful
    if( isRecoInter )
      if( fabs(higgs_inter[0].m() - 125.) < 40.0 && fabs(higgs_inter[1].m() - 125.) < 40.0 )
      {
        // Record cutflow
        FillHistogram("CF_inter", event_weight, 1.1);
        FillHistogram("CFN_inter", 1., 1.1);

        // Histograms for reconstructed Higgs candidates
        FillHistogram("pt_H0_inter_C1", event_weight, higgs_inter[0].pt());
        FillHistogram("pt_H1_inter_C1", event_weight, higgs_inter[1].pt());

        FillHistogram("m_H0_inter_C1", event_weight, higgs_inter[0].m());
        FillHistogram("m_H1_inter_C1", event_weight, higgs_inter[1].m());

        // Reconstruct di-Higgs system
        const fastjet::PseudoJet dihiggs_inter = higgs_inter[0] + higgs_inter[1];

        // Histograms for reconstructed di-Higgs system
        FillHistogram("m_HH_inter_C1", event_weight, dihiggs_inter.m());
        FillHistogram("pt_HH_inter_C1", event_weight, dihiggs_inter.pt());
        FillHistogram("dR_HH_inter_C1", event_weight, higgs_inter[0].delta_R(higgs_inter[1]) );
        FillHistogram("dPhi_HH_inter_C1", event_weight, getDPhi(higgs_inter[0].phi(), higgs_inter[1].phi()) );
        FillHistogram("dEta_HH_inter_C1", event_weight, fabs( higgs_inter[0].eta() - higgs_inter[1].eta()) );
      }
  }
  
  //=============================
  // C2: b-tagging
  //=============================
  if( nBBTaggedFatJets >= 2 )
  {
    // Higgs mass window cut
    if( fabs(bbFatJets[0].m() - 125.) < 40.0 && fabs(bbFatJets[1].m() - 125.) < 40.0 )
    {
      
      // b-tagging weights
      const double nB = nBSubjetsLR_vec[0] + nBSubjetsLR_vec[1];      // Number of true b-subjets
      const double nF = 4 - nB;  // Number of fake b-subjets

      // Reweighted event weight
      const double boost_weight = pow(btag_prob,nB)*pow(btag_mistag,nF)*event_weight;

    	// Record cutflow
    	FillHistogram("CF_boost", boost_weight, 2.1);
    	FillHistogram("CFN_boost", 1., 2.1);

    	// Histograms for reconstructed Higgs candidates
    	FillHistogram("pt_H0_boost_C2", boost_weight, bbFatJets[0].pt());
    	FillHistogram("pt_H1_boost_C2", boost_weight, bbFatJets[1].pt());
    	
    	FillHistogram("m_H0_boost_C2", boost_weight, bbFatJets[0].m());
    	FillHistogram("m_H1_boost_C2", boost_weight, bbFatJets[1].m());
    	
    	FillHistogram("ptHptH_boost_C2", boost_weight, bbFatJets[0].pt(), bbFatJets[1].pt());
    	FillHistogram("mHmH_boost_C2", boost_weight, bbFatJets[0].m(), bbFatJets[1].m());
    	
    	// Reconstruct di-Higgs system
    	const fastjet::PseudoJet dihiggs_boost = bbFatJets[0] + bbFatJets[1];
    	
    	// Histograms for reconstructed di-Higgs system
    	FillHistogram("m_HH_boost_C2", boost_weight, dihiggs_boost.m());
    	FillHistogram("pt_HH_boost_C2", boost_weight, dihiggs_boost.pt());
    	FillHistogram("dR_HH_boost_C2", boost_weight, bbFatJets[0].delta_R(bbFatJets[1]) );
    	FillHistogram("dPhi_HH_boost_C2", boost_weight, getDPhi(bbFatJets[0].phi(), bbFatJets[1].phi()) );
    	FillHistogram("dEta_HH_boost_C2", boost_weight, fabs( bbFatJets[0].eta() - bbFatJets[1].eta()) );
    	//std::cout << "dPhi " << getDPhi(bbFatJets[0].phi(), bbFatJets[1].phi()) << std::endl;
      }
  }
  if( nBTaggedJets >= 4 ){
    
    // Reconstruct Higgs candidates from b-tagged small-R jets
    std::vector<fastjet::PseudoJet> higgs_res;
    Reco_Resolved( bJets, higgs_res );
    
    // Higgs mass window cut
    if( fabs(higgs_res[0].m() - 125.) < 40.0 && fabs(higgs_res[1].m() - 125.) < 40.0 ){
        
      int nB = 0;
      for (int i=0; i<nBTaggedJets; i++)
        if (!isFakeSR_vec[i])
          nB++;

      const int nF = 4 - nB;


      // Reweighted event weight
      const double res_weight = pow(btag_prob,nB)*pow(btag_mistag,nF)*event_weight;

    	// Record cutflow
    	FillHistogram("CF_res", res_weight, 2.1);
    	FillHistogram("CFN_res", 1., 2.1);
          
    	// Histograms for reconstructed Higgs candidates
    	FillHistogram("pt_H0_res_C2", res_weight, higgs_res[0].pt());
    	FillHistogram("pt_H1_res_C2", res_weight, higgs_res[1].pt());
    	
    	FillHistogram("m_H0_res_C2", res_weight, higgs_res[0].m());
    	FillHistogram("m_H1_res_C2", res_weight, higgs_res[1].m());
    	
    	FillHistogram("ptHptH_res_C2", res_weight, higgs_res[0].pt(), higgs_res[1].pt());
    	FillHistogram("mHmH_res_C2", res_weight, higgs_res[0].m(), higgs_res[1].m());
    	
    	// Reconstruct di-Higgs system
    	const fastjet::PseudoJet dihiggs_res = higgs_res[0] + higgs_res[1];
    	
    	// Histograms for reconstructed di-Higgs system
    	FillHistogram("m_HH_res_C2", res_weight, dihiggs_res.m());
    	FillHistogram("pt_HH_res_C2", res_weight, dihiggs_res.pt());
    	FillHistogram("dR_HH_res_C2", res_weight, higgs_res[0].delta_R(higgs_res[1]) );
    	FillHistogram("dPhi_HH_res_C2", res_weight, getDPhi(higgs_res[0].phi(), higgs_res[1].phi()) );
    	FillHistogram("dEta_HH_res_C2", res_weight, fabs( higgs_res[0].eta() - higgs_res[1].eta()) );
    	//std::cout << "dPhi " << getDPhi(higgs_res[0].phi(), higgs_res[1].phi()) << std::endl;
    }
  }
  if( nBTaggedJets >= 2 &&  nBBTaggedFatJets == 1 )
  {
    // Reconstruct Higgs candidates from bb-tagged large-R and b-tagged small-R jets
    std::vector<fastjet::PseudoJet> higgs_inter;
    int nBJets_SR;
    bool isRecoInter = Reco_Intermediate( smallRJetsSel, isFakeSR_vec, largeRJetsSel[0], nBJets_SR, higgs_inter );
    
    // Check if reconstruction was successful
    if( isRecoInter )
      if( fabs(higgs_inter[0].m() - 125.) < 40.0 && fabs(higgs_inter[1].m() - 125.) < 40.0 )
      {
        const int nB = nBJets_SR + nBSubjetsLR_vec[0];
        const int nF = 4 - nB;

        if (nB > 4)
        {
          std::cerr << "ERROR: number of reconstructed b quarks > 4!"<<std::endl;
          exit(-1);
        }
        // Reweighted event weight
        const double inter_weight = pow(btag_prob,nB)*pow(btag_mistag,nF)*event_weight;

        // Record cutflow
        FillHistogram("CF_inter", inter_weight, 2.1);
        FillHistogram("CFN_inter", 1., 2.1);

        // Histograms for reconstructed Higgs candidates
        FillHistogram("pt_H0_inter_C2", inter_weight, higgs_inter[0].pt());
        FillHistogram("pt_H1_inter_C2", inter_weight, higgs_inter[1].pt());

        FillHistogram("m_H0_inter_C2", inter_weight, higgs_inter[0].m());
        FillHistogram("m_H1_inter_C2", inter_weight, higgs_inter[1].m());

        FillHistogram("ptHptH_inter_C2", inter_weight, higgs_inter[0].pt(), higgs_inter[1].pt());
        FillHistogram("mHmH_inter_C2", inter_weight, higgs_inter[0].m(), higgs_inter[1].m());

        // Reconstruct di-Higgs system
        const fastjet::PseudoJet dihiggs_inter = higgs_inter[0] + higgs_inter[1];

        // Histograms for reconstructed di-Higgs system
        FillHistogram("m_HH_inter_C2", inter_weight, dihiggs_inter.m());
        FillHistogram("pt_HH_inter_C2", inter_weight, dihiggs_inter.pt());
        FillHistogram("dR_HH_inter_C2", inter_weight, higgs_inter[0].delta_R(higgs_inter[1]) );
        FillHistogram("dPhi_HH_inter_C2", inter_weight, getDPhi(higgs_inter[0].phi(), higgs_inter[1].phi()) );
        FillHistogram("dEta_HH_inter_C2", inter_weight, fabs( higgs_inter[0].eta() - higgs_inter[1].eta()) );
        //std::cout << "dPhi " << getDPhi(higgs_inter[0].phi(), higgs_inter[1].phi()) << std::endl;
      }
  }
 
  // ************************************* MVA Output **********************************************************
  // Now save the ntuples to be used by the TMVA or the ANNs
  // In the UCL analysis they use
  //
  // m, y, pT of the 4b system and masses of the two dijets
  // 3 decay angles (in resp. rest frames) & 2 angles between decay planes
  // This is for the UCL-like strategy
  // sabe mass, pt and y of th 4b system
  // the two dijet masses
  // and all independent angular distances between the four b jets
  // totalNTuple<<"# signal source m4b pt4b y4b mHiggs1 mHiggs2 DeltaR_b1b2 DeltaR_b1b3 DeltaR_b1b4 DeltaR_b2b3 DeltaR_b2b4 DeltaR_b3b4 "<<std::endl;
  outputNTuple <<signal <<"\t"<<GetSample()<<"\t"<<event_weight<<std::endl;
  // Other combinations of kinematical variables could also be useful
  // Need to investigate the kinematics of the 4b final state in more detail
  // Pass event
  Pass(event_weight);

}


void OxfordCombinedRWAnalysis::BTagging( std::vector<fastjet::PseudoJet> const& jets_vec, std::vector<bool>& isFake_vec  )
{
  
  // Loop over all jets
  for( size_t i=0; i<jets_vec.size(); i++){

    // Get the jet constituents
    const std::vector<fastjet::PseudoJet>& jet_constituents = jets_vec[i].constituents();
    bool isbJet = false;

    // Loop over constituents and look for b quarks
    // b quarks must be above some minimum pt
    for(size_t j=0; j<jet_constituents.size(); j++)
    {
      // Flavour of jet constituent
      const int userid= jet_constituents.at(j).user_index();
      const double pt_bcandidate = jet_constituents.at(j).pt();

      if(abs(userid) == 5 && pt_bcandidate > pt_btagging)
        isbJet = true;    
    }

    // Push opposite of isBJet!
    isFake_vec.push_back(!isbJet);
  }
         
}

void OxfordCombinedRWAnalysis::BTaggingFJ( std::vector<fastjet::PseudoJet> const& largeRJets, std::vector<fastjet::PseudoJet> const& trackjets, std::vector<int>& nBSubJets_vec )
{
    // Loop over all fat jets
    for( size_t i=0; i<largeRJets.size(); i++)
    {
      // Get ghost associated track jets
      std::vector<fastjet::PseudoJet> subjets;
      get_assoc_trkjets( largeRJets.at(i), trackjets, subjets, false);

      const int nSubjets = std::min((int)subjets.size(), 2); // Restrict to leading 2 subjets
      int nBSubJets = 0;

      // Loop over subjets
      for( size_t j=0; j < nSubjets; j++)
      {
        bool isBJet = false;
        // Loop over constituents and look for b quarks
        // b quarks must be above some minimum pt 
        const std::vector<fastjet::PseudoJet>& subjet_constituents = sorted_by_pt( subjets[j].constituents() );
        for(size_t k=0; k < subjet_constituents.size(); k++)
        {
          // Flavour of jet constituent
          const int userid= subjet_constituents.at(k).user_index();
          const double pt_bcandidate = subjet_constituents.at(k).pt();

          if( abs(userid) == 5 && pt_bcandidate > pt_btagging )
            isBJet = true;
        }

        if (isBJet)
          nBSubJets++;
           
      }//end of loop over subjets
  	   
      nBSubJets_vec.push_back(nBSubJets);

    }//end of loop over jets
}


void OxfordCombinedRWAnalysis::Reco_Resolved( std::vector<fastjet::PseudoJet> const& bjets, std::vector<fastjet::PseudoJet>& higgs_vec ){

    // Get the pairing that minimizes |m_dj1 - m_dj2|
    double dijet_mass[4][4];
    for(int ijet=0;ijet<4;ijet++)
      for(int jjet=0;jjet<4;jjet++)
      {
        // Compute jet masses
        const fastjet::PseudoJet sum = bjets[ijet] + bjets[jjet];
        dijet_mass[ijet][jjet] = sum.m();
      }

    double mdj_diff_min = 1e20; // Some large number to begin
    int jet1_id1=10,jet1_id2=10,jet2_id1=10,jet2_id2=10;

    for(int ijet=0;ijet<4;ijet++)
      for(int jjet=ijet+1;jjet<4;jjet++)
        for(int ijet2=0;ijet2<4;ijet2++)
          for(int jjet2=ijet2+1;jjet2<4;jjet2++)
          {
            const double mdj1 = dijet_mass[ijet][jjet];
            const double mdj2 = dijet_mass[ijet2][jjet2];
            const double min_dj = fabs(mdj1 - mdj2);

            if(min_dj <  mdj_diff_min && ijet != ijet2  && jjet != jjet2 && jjet !=ijet2 && ijet != jjet2 )
            {
              mdj_diff_min = min_dj;
              jet1_id1 = ijet;
              jet1_id2 = jjet;
              jet2_id1 = ijet2;
              jet2_id2 = jjet2;
            }
          }

    // Construct the Higgs candidates
    fastjet::PseudoJet higgs1 = bjets.at( jet1_id1) + bjets.at( jet1_id2); 
    fastjet::PseudoJet higgs2 = bjets.at( jet2_id1) + bjets.at( jet2_id2);
    
    if( higgs1.pt() > higgs2.pt()){
      higgs_vec.push_back(higgs1);
      higgs_vec.push_back(higgs2);
    }
    else{
      higgs_vec.push_back(higgs2);
      higgs_vec.push_back(higgs1);    
    }
}


bool OxfordCombinedRWAnalysis::Reco_Intermediate( std::vector<fastjet::PseudoJet> const& bjets, 
                                                  std::vector<bool> const& isFakeSR_vec,
                                                  fastjet::PseudoJet const& fatjet, 
                                                  int& nBjets,
                                                  std::vector<fastjet::PseudoJet>& higgs_vec )
{
  
    // Identify small-R jets separated from merged Higgs
    std::vector<fastjet::PseudoJet> bjets_separated;
    std::vector<bool> bjets_separated_isFake;
    
    for( size_t i = 0; i < bjets.size(); i++ ){
      
      double dR = bjets.at(i).delta_R(fatjet);
      if( dR < 1.0 ) continue;
      
      bjets_separated.push_back( bjets.at(i) );
      bjets_separated_isFake.push_back( isFakeSR_vec[i] );
    }
    
    if( bjets_separated.size() < 2 ) return false;
    
    double mdj_diff_min = 1e20; // Some large number to begin
    int bjet_id1(100), bjet_id2(100);
    
    // Get the pairing that minimizes |m_dj1 - m_dj2|
    for(int ijet=0; ijet < (int)bjets_separated.size(); ijet++)
      for(int jjet=0; jjet < (int)bjets_separated.size(); jjet++)
      {
        // Compute jet masses
        const fastjet::PseudoJet sum = bjets_separated[ijet] + bjets_separated[jjet];
	
	double diff = fabs( sum.m() - fatjet.m() );
	if( diff < mdj_diff_min ){
	  
	  bjet_id1 = ijet;
	  bjet_id2 = jjet;
	  mdj_diff_min = diff;
	}
      }
    
    // Construct the Higgs candidates
    fastjet::PseudoJet higgs1 = fatjet; 
    fastjet::PseudoJet higgs2 = bjets_separated.at( bjet_id1) + bjets_separated.at( bjet_id2);
    nBjets = (bjets_separated_isFake[bjet_id1] ? 0:1) + (bjets_separated_isFake[bjet_id2] ? 0:1);
    
    if( higgs1.pt() > higgs2.pt()){
      higgs_vec.push_back(higgs1);
      higgs_vec.push_back(higgs2);
    }
    else{
      higgs_vec.push_back(higgs2);
      higgs_vec.push_back(higgs1);    
    }
    
    return true;
}