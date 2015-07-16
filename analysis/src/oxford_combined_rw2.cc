// oxford_res_fr.cc

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

#include "oxford_combined_rw2.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

#include <algorithm>

using namespace fastjet::contrib;

// pT cut for constituent b-quarks (GeV)
const double pt_btagging = 15;

// Boosted jet radius
const double BoostJetR=1.0;
const double GAjetR = 0.3; // Boosted subjet radius for ghost-association

// Resolved jet radius
const double ResJetR=0.4;

// beta-exponent for LST energy correlations
const double LST_beta = 2;

// Debugging
const bool debug = false;

// Exclusivity cut
const bool exclusive = true;

// Analysis settings
const int nAnalysis = 3;  const int nCuts = 7;
const std::string aString[nAnalysis] = {"_res", "_inter", "_boost"};
const std::string cString[nCuts] = {"_C0", "_C1a", "_C1b", "_C1c", "_C1d", "_C1e", "_C2"};


OxfordCombinedRW2Analysis::OxfordCombinedRW2Analysis(std::string const& sampleName):
Analysis("oxford_combined_rw2", sampleName)
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
  
  const double eta_min = -6.;
  const double eta_max = +6.;
  
  const double phi_min = -3.15;
  const double phi_max = +3.15;
  
  const double m_HH_min = 0.;
  const double m_HH_max = 600.; 
  
  const double pt_HH_min = 0.;
  const double pt_HH_max = 300.;
  
  const int nbins = 30;
  
  // ********************* Histogram definitions ******************

  for (int i=0; i< nAnalysis; i++)
  {
    BookHistogram(new YODA::Histo1D( nCuts, 0, nCuts ), "CF" + aString[i]);
    BookHistogram(new YODA::Histo1D( nCuts, 0, nCuts ), "CFN" + aString[i]);

    for (int j=0; j< nCuts; j++)
    {
      const std::string suffix = aString[i] + cString[j];
      
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_smallR" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_smallR" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_smallR" + suffix);
      
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_largeR" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_largeR" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_largeR" + suffix);
      
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1" + suffix);

      BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1" + suffix);
      
      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_H0" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_H1" + suffix);
      
      BookHistogram(new YODA::Histo1D(nbins, phi_min, phi_max), "phi_H0" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, phi_min, phi_max), "phi_H1" + suffix);

      BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH" + suffix);

      BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH" + suffix);

      BookHistogram(new YODA::Histo2D(nbins, pt_min, pt_max, nbins, pt_min, pt_max), "ptHptH" + suffix);
      BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mHmH" + suffix);
      
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_leadSJ_fj1" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_subleadSJ_fj1" + suffix);

      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_leadSJ_fj2" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_subleadSJ_fj2" + suffix);
      
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_leadSJ_fj" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_subleadSJ_fj" + suffix);

      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_leadSJ_fj1" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_subleadSJ_fj1" + suffix);

      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_leadSJ_fj2" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_subleadSJ_fj2" + suffix);
      
      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_leadSJ_fj" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_subleadSJ_fj" + suffix);
      
      BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj1" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj2" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj" + suffix);  

      BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj1" + suffix); 
      BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj2" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj" + suffix);      

      BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "C2_fj1" + suffix);  
      BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "C2_fj2" + suffix);
      BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "C2_fj" + suffix);   

      BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "D2_fj1" + suffix);  
      BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "D2_fj2" + suffix);
      BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "D2_fj" + suffix);       
      
      // Additional histograms from unreweighted analysis
      BookHistogram(new YODA::Histo1D(10, 0, 10), "N_SmallRJets" + suffix);
      BookHistogram(new YODA::Histo1D(10, 0, 10), "N_SmallRJets_BJets" + suffix);
      BookHistogram(new YODA::Histo1D(10, 0, 10), "N_SmallRJets_BTagged" + suffix);
    }
  }
    
  // ********************* Ntuple definition **********************

  const std::string tupleSpec = "# signal source weight pt_H0 pt_H1 pt_HH m_H0 m_H1 m_HH dR_HH dPhi_HH dEta_HH pt_H0_sub0 pt_H0_sub1 pt_H1_sub0 pt_H1_sub1";
  outputNTuple<<tupleSpec<<std::endl;

  const std::string root = "." + GetRoot() + GetSample() + "/";

  const std::string resDir = root+"resNTuple.dat";
  const std::string intDir = root+"intNTuple.dat";
  const std::string bstDir = root+"bstNTuple.dat";

  resNTuple.open(resDir.c_str());
  intNTuple.open(intDir.c_str());
  bstNTuple.open(bstDir.c_str());

  resNTuple << tupleSpec <<std::endl;
  intNTuple << tupleSpec <<" split12_fj tau21_fj C2_fj D2_fj"<<std::endl;
  bstNTuple << tupleSpec <<" split12_fj1 split12_fj2 tau21_fj1 tau21_fj2 C2_fj1 C2_fj2 D2_fj1 D2_fj2"<<std::endl;

  // Category overlap
  for (int i=0; i<nCuts; i++)
    BookHistogram(new YODA::Histo1D( 8, 0, 8 ), "Categories"+cString[i]);

}

void OxfordCombinedRW2Analysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
  Analysis::Analyse(signal, weightnorm, fs);

  // Set initial weight
  const double event_weight = weightnorm;
  double P_select_boost = 0; // Probability of boosted selection
  double P_select_inter = 0; // Probability of intermediate selection
  double P_select_resol = 0; // Probability of resolved selection

  // *************************************** Classification vectors *********************************

  std::vector<bool> bstClassified(nCuts,false);
  std::vector<bool> intClassified(nCuts,false);
  std::vector<bool> resClassified(nCuts,false);

  bstClassified[0] = true;
  intClassified[0] = true;
  resClassified[0] = true;

  // *************************************** General cut consts *************************************

  // Large-R jet cuts
  const double LR_minPT = 200;
  const double LR_maxEta = 2.0;

  // Small-R jet cuts
  const double SR_minPT = 40;
  const double SR_maxEta = 2.5;

  // Track jet cuts
  const double TJ_minPT = 50;
  const double TJ_maxEta = 2.5;

  // Higgs mass-window
  const double massWindow = 40;

  // ********************************* Jet clustering  ***************************************
  
  // Cluster small-R jets
  const fastjet::JetDefinition akt_res(fastjet::antikt_algorithm, ResJetR);
  const fastjet::ClusterSequence cs_akt_res(fs, akt_res);
  const std::vector<fastjet::PseudoJet> smallRJets_noCut = sorted_by_pt( cs_akt_res.inclusive_jets()  ); 
  const std::vector<fastjet::PseudoJet> smallRJets_pTcut = sorted_by_pt( cs_akt_res.inclusive_jets( SR_minPT )  ); 

  // Eta cut
  std::vector<fastjet::PseudoJet> smallRJets_etacut;
  for (size_t i=0; i<smallRJets_pTcut.size(); i++)
    if ( fabs(smallRJets_pTcut[i].eta()) <= SR_maxEta)
      smallRJets_etacut.push_back(smallRJets_pTcut[i]);

  // Cluster large-R jets
  const fastjet::JetDefinition akt_boost(fastjet::antikt_algorithm, BoostJetR);
  const fastjet::ClusterSequence cs_akt_bst(fs, akt_boost);
  const std::vector<fastjet::PseudoJet> largeRJets_noCut = sorted_by_pt( cs_akt_bst.inclusive_jets()  ); 
  const std::vector<fastjet::PseudoJet> largeRJets_pTcut = sorted_by_pt( cs_akt_bst.inclusive_jets( LR_minPT ) ); 

  // Eta cut
  std::vector<fastjet::PseudoJet> largeRJets_etacut;
  for (size_t i=0; i<largeRJets_pTcut.size(); i++)
  {
    if ( fabs(largeRJets_pTcut[i].eta()) <= LR_maxEta)
      largeRJets_etacut.push_back(largeRJets_pTcut[i]);
  }

  // Cluster small-R track jets 
  const fastjet::JetDefinition jd_subjets(fastjet::antikt_algorithm, GAjetR);
  const fastjet::ClusterSequence cs_subjets(fs, jd_subjets);
  const std::vector<fastjet::PseudoJet> trackjets_nocut = sorted_by_pt( cs_subjets.inclusive_jets()  );
  const std::vector<fastjet::PseudoJet> trackjets_pTcut = sorted_by_pt( cs_subjets.inclusive_jets( TJ_minPT ) );

  // Eta cut
  std::vector<fastjet::PseudoJet> trackJets_etacut;
  for (size_t i=0; i<trackjets_pTcut.size(); i++)
    if ( fabs(trackjets_pTcut[i].eta()) <= TJ_maxEta)
      trackJets_etacut.push_back(trackjets_pTcut[i]);

  // Final sorted jets
  const std::vector<fastjet::PseudoJet> smallRJets = sorted_by_pt(smallRJets_etacut);
  const std::vector<fastjet::PseudoJet> trackJets = sorted_by_pt(trackJets_etacut);

  // ********************************************* MDT *********************************************************

  // Check if jets are mass-drop tagged
  const fastjet::JetDefinition CA10(fastjet::cambridge_algorithm, 1.0);
  const fastjet::MassDropTagger md_tagger(mu, ycut);

  std::vector<fastjet::PseudoJet> MDTJets;
  for (size_t i=0; i<largeRJets_etacut.size(); i++)
  {
    const fastjet::ClusterSequence cs_sub( largeRJets_etacut[i].constituents(), CA10);
    const fastjet::PseudoJet ca_jet = sorted_by_pt(cs_sub.inclusive_jets())[0];
    const fastjet::PseudoJet tagged_jet = md_tagger(ca_jet);
    if ( tagged_jet != 0 )
      MDTJets.push_back(largeRJets_etacut[i]);
  }

  const std::vector<fastjet::PseudoJet> largeRJets = sorted_by_pt(MDTJets);

  // ***************************************** Initial histograms **********************************************

  FillHistogram("CF_res", event_weight, 0.1);
  FillHistogram("CF_inter", event_weight, 0.1);
  FillHistogram("CF_boost", event_weight, 0.1);

  FillHistogram("CFN_res", 1., 0.1);
  FillHistogram("CFN_inter", 1., 0.1);
  FillHistogram("CFN_boost", 1., 0.1);

  // ***************************************** B-Tagging **********************************************

  // b-tagging for large-R jets
  std::vector<int> nBSubjetsLR_vec;    // Vector specifying how many b subjets there are
  std::vector<int> nCSubjetsLR_vec;    // Vector specifying how many c subjets there are
  std::vector<int> nLSubjetsLR_vec;    // Vector specifying how many l subjets there are
  std::vector<fastjet::PseudoJet> leading_subjet;
  std::vector<fastjet::PseudoJet> subleading_subjet;

  BTagging( largeRJets, trackJets, leading_subjet, subleading_subjet, nBSubjetsLR_vec, nCSubjetsLR_vec, nLSubjetsLR_vec);
  if( largeRJets.size() != nBSubjetsLR_vec.size() )
    std::cout << "ERROR: b-tagging vector sizes don't match number of fat jets" << std::endl;

  // b-tagging for small-R jets
  std::vector<btagType> tagType_SR;
  BTagging( smallRJets, tagType_SR );

  // **************************************** Boosted analysis *********************************************

  // Boosted initial histograms
  if (largeRJets_noCut.size() >= 2) // Clustering
  {
    FillHistogram("CF_boost", event_weight, 1.1);
    FillHistogram("CFN_boost", 1., 1.1);
    JetFill(  smallRJets, largeRJets_noCut, "boost", 1, event_weight );
    bstClassified[1] = true;

    if (largeRJets_pTcut.size() >= 2) // pT cut
    {
      FillHistogram("CF_boost", event_weight, 2.1);
      FillHistogram("CFN_boost", 1., 2.1);
      JetFill(  smallRJets, largeRJets_pTcut, "boost", 2, event_weight );
      bstClassified[2] = true;

      if (largeRJets_etacut.size() >= 2) // Eta cut
      {
        FillHistogram("CF_boost", event_weight, 3.1);
        FillHistogram("CFN_boost", 1., 3.1);
        JetFill(  smallRJets, largeRJets_etacut, "boost", 3, event_weight );
        bstClassified[3] = true;

        if( largeRJets.size() >= 2 ) // MDT
        {
          HiggsFill(largeRJets[0], largeRJets[1], "boost", 4, event_weight);
          JetFill(  smallRJets, largeRJets, "boost", 4, event_weight );
          bstClassified[4] = true;

          // Higgs mass-window
          const double diffHiggs_0 = fabs(largeRJets[0].m() - 125.);
          const double diffHiggs_1 = fabs(largeRJets[1].m() - 125.);

          if( (diffHiggs_0 < massWindow) && (diffHiggs_1 < massWindow) )
          {
            HiggsFill(largeRJets[0], largeRJets[1], "boost", 5, event_weight);
            BoostFill(largeRJets[0], largeRJets[1], "boost", 5, event_weight);
            JetFill(  smallRJets, largeRJets, "boost", 5, event_weight );
            SubJetFill( leading_subjet, subleading_subjet, "boost", 5, event_weight);
            bstClassified[5] = true;

            // b-tagging weights
            const double nB = nBSubjetsLR_vec[0] + nBSubjetsLR_vec[1];      // Number of true b-subjets
            const double nC = nCSubjetsLR_vec[0] + nCSubjetsLR_vec[1];      // Number of fake b-subjets
            const double nL = nLSubjetsLR_vec[0] + nLSubjetsLR_vec[1];      // Number of fake b-subjets

            // Selection probability
            if (nB+nC+nL == 4) // Selected 4 candidates
            {
              P_select_boost = btagProb(4,nB,nC,nL);

              // Reweighted event weight
              const double boost_weight = P_select_boost*event_weight;
              const fastjet::PseudoJet dihiggs_boost = largeRJets[0] + largeRJets[1];

              HiggsFill(largeRJets[0], largeRJets[1], "boost", 6, boost_weight);
              BoostFill(largeRJets[0], largeRJets[1], "boost", 6, boost_weight);
              JetFill(  smallRJets, largeRJets, "boost", 6, boost_weight );
              SubJetFill( leading_subjet, subleading_subjet, "boost", 6, boost_weight);
              bstClassified[6] = true;

              // Calculate some substructure variables
              const std::vector<double> split12_vec = SplittingScales( largeRJets );
              const std::vector<double> tau21_vec = NSubjettiness( largeRJets, BoostJetR );

              // C2 energy correlation double-ratio
              const double C2_fj1 = LST_C2(LST_beta, largeRJets[0]);
              const double C2_fj2 = LST_C2(LST_beta, largeRJets[1]);

              // D2 energy correlation double-ratio
              const double D2_fj1 = LMN_D2(LST_beta, largeRJets[0]);
              const double D2_fj2 = LMN_D2(LST_beta, largeRJets[1]);

              // Fill tuple
              bstNTuple << signal <<"\t"<<GetSample()<<"\t"<<boost_weight << "\t"
                  << largeRJets[0].pt() << "\t"
                  << largeRJets[1].pt() << "\t"
                  << dihiggs_boost.pt() << "\t"
                  << largeRJets[0].m() << "\t"
                  << largeRJets[1].m() << "\t"
                  << dihiggs_boost.m() << "\t"
                  << largeRJets[0].delta_R(largeRJets[1])  << "\t"
                  << getDPhi(largeRJets[0].phi(), largeRJets[1].phi())  << "\t"
                  << fabs( largeRJets[0].eta() - largeRJets[1].eta())  << "\t"
                  << leading_subjet[0].pt() << "\t"
                  << subleading_subjet[0].pt() << "\t"
                  << leading_subjet[1].pt() << "\t"
                  << subleading_subjet[1].pt() << "\t"
                  << split12_vec[0] << "\t"
                  << split12_vec[1] << "\t"
                  << tau21_vec[0] << "\t"
                  << tau21_vec[1] << "\t"
                  << C2_fj1 << "\t"
                  << C2_fj2 << "\t"
                  << D2_fj1 << "\t"
                  << D2_fj2 << "\t"
                  <<std::endl;

              Pass(boost_weight);
              Cut("BoostedCut", event_weight - boost_weight );
            }
          }
        }
      }
    }
  }

  // ************************************* Intermediate analysis ********************************************
/*
  if (smallRJets_noCut.size() >= 2 &&  largeRJets_noCut.size() == 1 ) // Clustering
  {
    FillHistogram("CF_inter", event_weight, 1.1);
    FillHistogram("CFN_inter", 1., 1.1);
    JetFill(  smallRJets_noCut, largeRJets_noCut, "inter", 1, event_weight );
    intClassified[1] = true;
  }

  if (smallRJets_pTcut.size() >= 2 && largeRJets_pTcut.size() == 1) // pT cut
  {
    FillHistogram("CF_inter", event_weight, 2.1);
    FillHistogram("CFN_inter", 1., 2.1);
    JetFill(  smallRJets_pTcut, largeRJets_pTcut, "inter", 2, event_weight );
    intClassified[2] = true;
  }

  if (smallRJets_etacut.size() >= 2 && largeRJets_etacut.size() == 1) // Eta cut
  {
    FillHistogram("CF_inter", event_weight, 3.1);
    FillHistogram("CFN_inter", 1., 3.1);
    JetFill(  smallRJets_etacut, largeRJets_etacut, "inter", 3, event_weight );
    intClassified[3] = true;
  }

  if( smallRJets.size() >= 2 &&  largeRJets.size() == 1 ) // MDT + reco cut
  {
    // Reconstruct Higgs candidates from large-R and small-R jets
    std::vector<fastjet::PseudoJet> higgs_inter; int nBJets_SR = 0; int nFJets_SR = 0;
    fastjet::PseudoJet res_leading_subjet; fastjet::PseudoJet res_subleading_subjet; 
    const bool isRecoInter = Reco_Intermediate( smallRJets, isTrueSR_vec, isFakeSR_vec, largeRJets[0], nBJets_SR, nFJets_SR, res_leading_subjet, res_subleading_subjet, higgs_inter );

    if( isRecoInter )
    {   
      HiggsFill(higgs_inter[0], higgs_inter[1], "inter", 4, event_weight);
      BoostFill(largeRJets[0], "inter", 4, event_weight);
      JetFill(  smallRJets, largeRJets, "inter", 4, event_weight );
      intClassified[4] = true;

      // Higgs mass-window
      const double diffHiggs_0 = fabs(higgs_inter[0].m() - 125.);
      const double diffHiggs_1 = fabs(higgs_inter[1].m() - 125.);

      if( (diffHiggs_0 < massWindow) && (diffHiggs_1 < massWindow) )
      {
        HiggsFill(higgs_inter[0], higgs_inter[1], "inter", 5, event_weight);
        BoostFill(largeRJets[0], "inter", 5, event_weight);
        JetFill(  smallRJets, largeRJets, "inter", 5, event_weight );
        intClassified[5] = true;

        // Determine number of fake bJets
        const int nB = nBJets_SR + nBSubjetsLR_vec[0];
        const int nF = nFJets_SR + nFSubjetsLR_vec[0];

        if (nB+nF == 4)
        {
          // Selection probability
          P_select *= pow(btag_prob,nB)*pow(btag_mistag,nF);

          // Reweighted event weight
          const double inter_weight = pow(btag_prob,nB)*pow(btag_mistag,nF)*event_weight;
          const fastjet::PseudoJet dihiggs_inter = higgs_inter[0] + higgs_inter[1];

          HiggsFill(higgs_inter[0], higgs_inter[1], "inter", 6, inter_weight);
          BoostFill(largeRJets[0], "inter", 6, inter_weight);
          JetFill(  smallRJets, largeRJets, "inter", 6, inter_weight );
          SubJetFill( leading_subjet, subleading_subjet, "inter", 6, inter_weight);
          intClassified[6] = true;

          if (!selected)
          {
            // Exclusivity cut
            HiggsFill(higgs_inter[0], higgs_inter[1], "inter", 7, inter_weight);
            BoostFill(largeRJets[0], "inter", 7, inter_weight);
            JetFill(  smallRJets, largeRJets, "inter", 7, inter_weight );
            SubJetFill( leading_subjet, subleading_subjet, "inter", 7, inter_weight);
            intClassified[7] = true;

            // Calculate some substructure variables
            const double split12 = SplittingScales( largeRJets[0] );
            const double tau21 = NSubjettiness( largeRJets[0], BoostJetR );
            const double C2 = LST_C2(LST_beta, largeRJets[0]);
            const double D2 = LMN_D2(LST_beta, largeRJets[0]);

            // Fill tuple
            intNTuple << signal <<"\t"<<GetSample()<<"\t"<<inter_weight << "\t"
                      << higgs_inter[0].pt() << "\t"
                      << higgs_inter[1].pt() << "\t"
                      << dihiggs_inter.pt() << "\t"
                      << higgs_inter[0].m() << "\t"
                      << higgs_inter[1].m() << "\t"
                      << dihiggs_inter.m() << "\t"
                      << higgs_inter[0].delta_R(higgs_inter[1]) << "\t"
                      << getDPhi(higgs_inter[0].phi(), higgs_inter[1].phi()) << "\t"
                      << fabs( higgs_inter[0].eta() - higgs_inter[1].eta())  << "\t"
                      << leading_subjet[0].pt() << "\t"
                      << subleading_subjet[0].pt() << "\t"
                      << res_leading_subjet.pt() << "\t"
                      << res_subleading_subjet.pt() << "\t"
                      << split12 << "\t"
                      << tau21 << "\t"
                      << C2 << "\t"
                      << D2 << "\t"
                      <<std::endl;

            // Final
            Pass(inter_weight);  selected = true;
            Cut("IntermediateCut", event_weight - inter_weight );
          }
        }
      }
    }
  }
*/

  // ************************************* Resolved analysis ********************************************

  if (smallRJets_noCut.size() >= 4) // Clustering
  {
    FillHistogram("CF_res", event_weight, 1.1);
    FillHistogram("CFN_res", 1., 1.1);
    JetFill(  smallRJets_noCut, largeRJets, "res", 1, event_weight );
    resClassified[1] = true;

    if (smallRJets_pTcut.size() >= 4) // pT cut
    {
      FillHistogram("CF_res", event_weight, 2.1);
      FillHistogram("CFN_res", 1., 2.1);
      JetFill(  smallRJets_pTcut, largeRJets, "res", 2, event_weight );
      resClassified[2] = true;

      if (smallRJets_etacut.size() >= 4) // Eta cut
      {
        FillHistogram("CF_res", event_weight, 3.1);
        FillHistogram("CFN_res", 1., 3.1);
	      JetFill(  smallRJets_etacut, largeRJets, "res", 3, event_weight );
        resClassified[3] = true;

        if( smallRJets.size() >= 4 )
        {
          // Reconstruct Higgs candidates from small-R jets
          std::vector<fastjet::PseudoJet> higgs_res;
          std::vector<fastjet::PseudoJet> higgs0_res;
          std::vector<fastjet::PseudoJet> higgs1_res;

          Reco_Resolved( smallRJets, higgs_res, higgs0_res, higgs1_res );
          HiggsFill( higgs_res[0], higgs_res[1], "res", 4, event_weight );
	        JetFill(  smallRJets, largeRJets, "res", 4, event_weight );
          resClassified[4] = true;

          // Higgs mass-window
          const double diffHiggs_0 = fabs(higgs_res[0].m() - 125.);
          const double diffHiggs_1 = fabs(higgs_res[1].m() - 125.);

          if( ( diffHiggs_0 < massWindow ) && ( diffHiggs_1 < massWindow ) )
          {
            HiggsFill( higgs_res[0], higgs_res[1], "res", 5, event_weight );
	          JetFill(  smallRJets, largeRJets, "res", 5, event_weight );
            resClassified[5] = true;

            // Determine number of real and fake b-jets
            int nB = 0; int nC = 0; int nL = 0;
            for (int i=0; i<4; i++)
              switch (tagType_SR[i])
              {
                case NTAG:
                  break;

                case BTAG:
                  nB++;
                  break;   

                case CTAG:
                  nC++;
                  break; 

                case LTAG:
                  nL++;
                  break; 
              }

            if (nB+nC+nL == 4) // Need 4 b-tags
            {
              // Selection probability
              P_select_resol = btagProb( 4, nB, nC, nL);

              // Reweighted event weight
              const double P_exclusive = exclusive ? (1.0-P_select_inter)*(1.0-P_select_boost):1;
              const double res_weight = P_select_resol*P_exclusive*event_weight;
              const fastjet::PseudoJet dihiggs_res = higgs_res[0] + higgs_res[1];

              HiggsFill(higgs_res[0], higgs_res[1], "res", 6, res_weight); 
              JetFill(  smallRJets, largeRJets, "res", 6, res_weight );
              resClassified[6] = true;

              resNTuple << signal <<"\t"<<GetSample()<<"\t"<<res_weight << "\t"
                        << higgs_res[0].pt() << "\t"
                        << higgs_res[1].pt() << "\t"
                        << dihiggs_res.pt() << "\t"
                        << higgs_res[0].m() << "\t"
                        << higgs_res[1].m() << "\t"
                        << dihiggs_res.m() << "\t"
                        << higgs_res[0].delta_R(higgs_res[1]) << "\t"
                        << getDPhi(higgs_res[0].phi(), higgs_res[1].phi()) << "\t"
                        << fabs( higgs_res[0].eta() - higgs_res[1].eta())  << "\t"
                        << higgs0_res[0].pt() << "\t"
                        << higgs0_res[1].pt() << "\t"
                        << higgs1_res[0].pt() << "\t"
                        << higgs1_res[1].pt() << "\t"
                        <<std::endl;

              Pass(res_weight); 
              Cut("ResolvedCut", event_weight - res_weight);
            }
          }
        }
      }
    }
  }

  // ************************************* Categorisation ********************************************
 
  for (int icut=0; icut<nCuts; icut++)
  {
    const std::string histoname = "Categories"+cString[icut];
    double coord = 0; double P = 1;


    if (resClassified[icut] && !intClassified[icut] && !bstClassified[icut])      { coord = 0.1; P=P_select_resol; }
    else if (!resClassified[icut] && intClassified[icut] && !bstClassified[icut]) { coord = 1.1; P=P_select_inter; }
    else if (!resClassified[icut] && !intClassified[icut] && bstClassified[icut]) { coord = 2.1; P=P_select_boost; }
    else if (resClassified[icut] && intClassified[icut] && !bstClassified[icut])  { coord = 3.1; P=P_select_resol*P_select_inter; }
    else if (resClassified[icut] && !intClassified[icut] && bstClassified[icut])  { coord = 4.1; P=P_select_resol*P_select_boost; }
    else if (!resClassified[icut] && intClassified[icut] && bstClassified[icut])  { coord = 5.1; P=P_select_inter*P_select_boost; }
    else if (resClassified[icut] && intClassified[icut] && bstClassified[icut])   { coord = 6.1; P=P_select_resol*P_select_inter*P_select_boost; }

    FillHistogram(histoname, (icut < 6) ? event_weight:P*event_weight, coord);
  }

  return Cut ("Uncategorised", (1.0-P_select_resol)*(1.0-P_select_inter)*(1.0-P_select_boost)*event_weight);
}

// Small-R B-tagging
void OxfordCombinedRW2Analysis::BTagging( std::vector<fastjet::PseudoJet> const& jets_vec, std::vector<btagType>& btag_vec  )
{
  // Loop over all jets
  for( size_t i=0; i<jets_vec.size(); i++)
  {
    // Initial categorisation
    btagType type = NTAG;

    // Loop over constituents and look for b quarks
    const std::vector<fastjet::PseudoJet>& jet_constituents = jets_vec[i].constituents();
    for(size_t j=0; j<jet_constituents.size(); j++)
    {
      // Flavour of jet constituent
      const int userid= jet_constituents.at(j).user_index();
      const double pt_bcandidate = jet_constituents.at(j).pt();

      if (pt_bcandidate > pt_btagging)
      {
        if(abs(userid) == 5) type = BTAG;
        if(abs(userid) == 4 && type != BTAG ) type = CTAG;
        if( type == NTAG ) type = LTAG;
      }
    }

    btag_vec.push_back(type);
  }
         
}

void OxfordCombinedRW2Analysis::BTagging( std::vector<fastjet::PseudoJet> const& largeRJets,
                                          std::vector<fastjet::PseudoJet> const& trackjets,
                                          std::vector<fastjet::PseudoJet>& subjets1,
                                          std::vector<fastjet::PseudoJet>& subjets2,
                                          std::vector<int>& nBSubJets_vec,
                                          std::vector<int>& nCSubJets_vec,
                                          std::vector<int>& nLSubJets_vec  )
{
    // Loop over all fat jets
    for( size_t i=0; i<largeRJets.size(); i++)
    {
      // Get ghost associated track jets
      std::vector<fastjet::PseudoJet> subjets;
      get_assoc_trkjets( largeRJets.at(i), trackjets, subjets, false);

      // Failed to cluster 2 subjets
      subjets = sorted_by_pt(subjets);
      if (subjets.size() < 2)
      {
        subjets1.push_back(fastjet::PseudoJet());
        subjets2.push_back(fastjet::PseudoJet());
      }
      else
      {
        subjets1.push_back(subjets[0]);
        subjets2.push_back(subjets[1]);
      }

      // B-tag subjets
      std::vector<btagType> btag_vec;
      BTagging( subjets, btag_vec );

      const int nSub = std::min((int)subjets.size(), 2); // Restrict to leading 2 subjets
      const int nBSubJets = std::count(btag_vec.begin(), btag_vec.begin() + nSub, BTAG);  // Number of b's
      const int nCSubJets = std::count(btag_vec.begin(), btag_vec.begin() + nSub, CTAG);  // Number of c's
      const int nLSubJets = std::count(btag_vec.begin(), btag_vec.begin() + nSub, LTAG);  // Number of u,d,s,g
      	   
      nBSubJets_vec.push_back(nBSubJets);
      nCSubJets_vec.push_back(nCSubJets);
      nLSubJets_vec.push_back(nLSubJets);

    }//end of loop over jets
}


void OxfordCombinedRW2Analysis::Reco_Resolved( std::vector<fastjet::PseudoJet> const& bjets, // Input b-jets
                                              std::vector<fastjet::PseudoJet>& higgs_vec,   // Returned Higgs candidates
                                              std::vector<fastjet::PseudoJet>& higgs0_vec,  // Leading higgs subjets
                                              std::vector<fastjet::PseudoJet>& higgs1_vec ) // Subleading higgs subjets
{

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
    
    if( higgs1.pt() > higgs2.pt())
    {
      higgs_vec.push_back(higgs1);
      higgs_vec.push_back(higgs2);

      higgs0_vec.push_back(bjets.at( jet1_id1));
      higgs0_vec.push_back(bjets.at( jet1_id2));

      higgs1_vec.push_back(bjets.at( jet2_id1));
      higgs1_vec.push_back(bjets.at( jet2_id2));

    }
    else{
      higgs_vec.push_back(higgs2);
      higgs_vec.push_back(higgs1);  

      higgs1_vec.push_back(bjets.at( jet1_id1));
      higgs1_vec.push_back(bjets.at( jet1_id2));

      higgs0_vec.push_back(bjets.at( jet2_id1));
      higgs0_vec.push_back(bjets.at( jet2_id2));  
    }

    // Sort vectors
    higgs0_vec = sorted_by_pt(higgs0_vec);
    higgs1_vec = sorted_by_pt(higgs1_vec);

}


bool OxfordCombinedRW2Analysis::Reco_Intermediate( std::vector<fastjet::PseudoJet> const& bjets, 
                                                  std::vector<bool> const& isTrueSR_vec,
                                                  std::vector<bool> const& isFakeSR_vec,
                                                  fastjet::PseudoJet const& fatjet, 
                                                  int& nBjets,
                                                  int& nFjets,
                                                  fastjet::PseudoJet& lead_subjet,
                                                  fastjet::PseudoJet& sublead_subjet,
                                                  std::vector<fastjet::PseudoJet>& higgs_vec )
{
  // Identify small-R jets separated from merged Higgs
  std::vector<fastjet::PseudoJet> bjets_separated;
  std::vector<bool> bjets_separated_isFake;
  std::vector<bool> bjets_separated_isTrue;

  for( size_t i = 0; i < bjets.size(); i++ )
  {  
    double dR = bjets.at(i).delta_R(fatjet);
    if( dR < 1.2 ) continue;
    
    bjets_separated.push_back( bjets.at(i) );
    bjets_separated_isFake.push_back( isFakeSR_vec[i] );
    bjets_separated_isTrue.push_back( isTrueSR_vec[i] );
  }
  
  if( bjets_separated.size() < 2 ) return false;
  
  double mdj_diff_min = 1e20; // Some large number to begin
  int bjet_id1(100), bjet_id2(100);
  
  int n_sep_jet_candidates = 2;
//   int n_sep_jet_candidates = (int)bjets_separated.size()
  
  // Get the pairing that minimizes |m_dj1 - m_dj2|
  for(int ijet=0; ijet < n_sep_jet_candidates; ijet++)
    for(int jjet=0; jjet < n_sep_jet_candidates; jjet++)
    {
      // Compute jet masses
      const fastjet::PseudoJet sum = bjets_separated[ijet] + bjets_separated[jjet];

      double diff = fabs( sum.m() - fatjet.m() );
      if( diff < mdj_diff_min )
      {
        bjet_id1 = ijet;
        bjet_id2 = jjet;
        mdj_diff_min = diff;
      }
    }
  
  // Construct the Higgs candidates
  fastjet::PseudoJet higgs1 = fatjet; 
  fastjet::PseudoJet higgs2 = bjets_separated.at( bjet_id1) + bjets_separated.at( bjet_id2);
  nBjets = (bjets_separated_isTrue[bjet_id1] ? 1:0) + (bjets_separated_isTrue[bjet_id2] ? 1:0);
  nFjets = (bjets_separated_isFake[bjet_id1] ? 1:0) + (bjets_separated_isFake[bjet_id2] ? 1:0);

  if( higgs1.pt() > higgs2.pt()){
    higgs_vec.push_back(higgs1);
    higgs_vec.push_back(higgs2);
  }
  else{
    higgs_vec.push_back(higgs2);
    higgs_vec.push_back(higgs1);    
  }

  // Subjets
  if (bjets_separated.at(bjet_id1).pt() > bjets_separated.at(bjet_id2).pt())
  {
    lead_subjet = bjets_separated.at( bjet_id1);
    sublead_subjet = bjets_separated.at( bjet_id2);
  }
  else
  {
    sublead_subjet = bjets_separated.at( bjet_id1);
    lead_subjet = bjets_separated.at( bjet_id2);
  }
  
  return true;
}

// Fill basic jet quantities
void OxfordCombinedRW2Analysis::JetFill(  std::vector<fastjet::PseudoJet> const& smallRJets,
					   std::vector<fastjet::PseudoJet> const& largeRJets,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight )
{
  
    // Histo fill suffix
    const std::string suffix = "_" + analysis + cString[cut];

    // Fill small-R jets
    for( int i = 0; i < (int) smallRJets.size(); i++ ){
      
	FillHistogram("pt_smallR" + suffix, weight, smallRJets[i].pt());
	FillHistogram("eta_smallR" + suffix, weight, smallRJets[i].eta());
	FillHistogram("m_smallR" + suffix, weight, smallRJets[i].m());
    }

    // Fill large-R jets
    for( int i = 0; i < (int) largeRJets.size(); i++ ){
      
	FillHistogram("pt_largeR" + suffix, weight, largeRJets[i].pt());
	FillHistogram("eta_largeR" + suffix, weight, largeRJets[i].eta());
	FillHistogram("m_largeR" + suffix, weight, largeRJets[i].m());
    }
}


// Fill basic subjet quantities
void OxfordCombinedRW2Analysis::SubJetFill(  std::vector<fastjet::PseudoJet> const& leading_subjet,
					      std::vector<fastjet::PseudoJet> const& subleading_subjet,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight )
{
  
    // Histo fill suffix
    const std::string suffix = "_" + analysis + cString[cut];
  
    // Fill histograms for boosted
    if ( analysis == "boost" ){
      
	if ( leading_subjet.size() < 2 )
	    std::cerr << "SubJetFill WARNING: Less than two leading subjets "<<analysis<<"  "<<cut<<"  "<< leading_subjet.size() <<std::endl;
	if ( subleading_subjet.size() < 2 )
	    std::cerr << "SubJetFill WARNING: Less than two subleading subjets "<<analysis<<"  "<<cut<<"  "<< subleading_subjet.size() <<std::endl;
	
	FillHistogram("pt_leadSJ_fj1" + suffix, weight, leading_subjet[0].pt());
	FillHistogram("pt_subleadSJ_fj1" + suffix, weight, subleading_subjet[0].pt());
	
	FillHistogram("pt_leadSJ_fj2" + suffix, weight, leading_subjet[1].pt());
	FillHistogram("pt_subleadSJ_fj2" + suffix, weight, subleading_subjet[1].pt());

	FillHistogram("eta_leadSJ_fj1" + suffix, weight, leading_subjet[0].eta());
	FillHistogram("eta_subleadSJ_fj1" + suffix, weight, subleading_subjet[0].eta());
	
	FillHistogram("eta_leadSJ_fj2" + suffix, weight, leading_subjet[1].eta());
	FillHistogram("eta_subleadSJ_fj2" + suffix, weight, subleading_subjet[1].eta());
	
    }
    
    // Fill histograms for intermediate
    if ( analysis == "inter" ){
      
	if ( leading_subjet.size() < 1 )
	    std::cerr << "SubJetFill WARNING: Less than one leading subjet"<<analysis<<"  "<<cut<<"  "<< leading_subjet.size() <<std::endl;
	if ( subleading_subjet.size() < 1 )
	    std::cerr << "SubJetFill WARNING: Less than one subleading subjet "<<analysis<<"  "<<cut<<"  "<< subleading_subjet.size() <<std::endl;
	
	FillHistogram("pt_leadSJ_fj" + suffix, weight, leading_subjet[0].pt());
	FillHistogram("pt_subleadSJ_fj" + suffix, weight, subleading_subjet[0].pt());
	
	FillHistogram("eta_leadSJ_fj" + suffix, weight, leading_subjet[0].eta());
	FillHistogram("eta_subleadSJ_fj" + suffix, weight, subleading_subjet[0].eta());
    }
    
}

// General fill for reconstructed higgs quantities
void OxfordCombinedRW2Analysis::HiggsFill(fastjet::PseudoJet const& H0,
                                         fastjet::PseudoJet const& H1,
                                         std::string const& analysis, 
                                         size_t const& cut, 
                                         double const& weight)
{
  if (H0.pt() < H1.pt())
    std::cerr << "HiggsFill WARNING: pT ordering incorrect! "<<analysis<<"  "<<cut<<"  "<<H0.pt() << "  "<<H1.pt()<<std::endl;

  if( debug ) std::cout << "HiggsFill INFO: cut = " << cut << std::endl;
  if( debug ) std::cout << "HiggsFill INFO: analysis = " << analysis << std::endl;

  // Histo fill suffix
  const std::string suffix = "_" + analysis + cString[cut];
    
  // Record cutflow
  FillHistogram("CF_" +analysis, weight, cut + 0.1);
  FillHistogram("CFN_"+analysis, 1., cut + 0.1);

  // Histograms for reconstructed Higgs candidates
  FillHistogram("pt_H0" + suffix, weight, H0.pt());
  FillHistogram("pt_H1" + suffix, weight, H1.pt());

  FillHistogram("m_H0" + suffix, weight, H0.m());
  FillHistogram("m_H1" + suffix, weight, H1.m());
  
  FillHistogram("eta_H0" + suffix, weight, H0.eta());
  FillHistogram("eta_H1" + suffix, weight, H1.eta());
  
  FillHistogram("phi_H0" + suffix, weight, H0.phi());
  FillHistogram("phi_H1" + suffix, weight, H1.phi());

  FillHistogram("ptHptH" + suffix, weight, H0.pt(), H1.pt());
  FillHistogram("mHmH" + suffix, weight, H0.m(), H1.m());

  FillHistogram("dR_HH" + suffix, weight, H0.delta_R(H1) );
  FillHistogram("dPhi_HH" + suffix, weight, getDPhi(H0.phi(), H1.phi()) );
  FillHistogram("dEta_HH" + suffix, weight, fabs( H0.eta() - H1.eta()) );

  // Reconstruct di-Higgs system
  const fastjet::PseudoJet dihiggs = H0 + H1;

  // Histograms for reconstructed di-Higgs system
  FillHistogram("m_HH" + suffix, weight, dihiggs.m());
  FillHistogram("pt_HH" + suffix, weight, dihiggs.pt());

}


void OxfordCombinedRW2Analysis::BoostFill( fastjet::PseudoJet const& H,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight )
{
  // Histo fill suffix
  const std::string suffix = "_" + analysis + cString[cut];

  // Splitting scales
  const double split12 = SplittingScales( H );

  // 2-subjettiness / 1-subjettiness
  const double tau21 = NSubjettiness( H, BoostJetR );

  // C2/D2 energy correlation double-ratio
  const double C2 = LST_C2(LST_beta, H);
  const double D2 = LMN_D2(LST_beta, H);

  FillHistogram("split12_fj" + suffix, weight, split12);
  FillHistogram("tau21_fj" + suffix, weight, tau21);
  FillHistogram("C2_fj" + suffix, weight, C2);
  FillHistogram("D2_fj" + suffix, weight, D2);
}

void OxfordCombinedRW2Analysis::BoostFill( fastjet::PseudoJet const& H0,
                                          fastjet::PseudoJet const& H1,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight )
{
  if (H0.pt() < H1.pt())
    std::cerr << "BoostFill WARNING: pT ordering incorrect! "<<analysis<<"  "<<cut<<"  "<<H0.pt() << "  "<<H1.pt()<<std::endl;
  
  if( debug ) std::cout << "BoostFill INFO: cut = " << cut << std::endl;
  if( debug ) std::cout << "BoostFill INFO: analysis = " << analysis << std::endl;

  // Histo fill suffix
  const std::string suffix = "_" + analysis + cString[cut];

  // Splitting scales
  if( debug ) std::cout << "BoostFill INFO: Calculate splitting scales" << std::endl;
  if( debug ) std::cout << "BoostFill INFO! H0 pt = " << H0.pt() << " H1 pt = "<< H1.pt() << std::endl;
  const double split12_fj1 = SplittingScales( H0 );
  const double split12_fj2 = SplittingScales( H1 );

  // 2-subjettiness / 1-subjettiness
  if( debug ) std::cout << "BoostFill INFO: Calculate n-subjettiness" << std::endl;
  const double tau21_fj1 = NSubjettiness( H0, BoostJetR );
  const double tau21_fj2 = NSubjettiness( H1, BoostJetR );

  // C2 energy correlation double-ratio
  if( debug ) std::cout << "BoostFill INFO: Calculate C2" << std::endl;
  const double C2_fj1 = LST_C2(LST_beta, H0);
  const double C2_fj2 = LST_C2(LST_beta, H1);

  // D2 energy correlation double-ratio
  if( debug ) std::cout << "BoostFill INFO: Calculate D2" << std::endl;
  const double D2_fj1 = LMN_D2(LST_beta, H0);
  const double D2_fj2 = LMN_D2(LST_beta, H1);
  
  if( debug ) std::cout << "BoostFill INFO! H0 split12 = " << split12_fj1 << " H1 split12 = "<< split12_fj2 << std::endl;
  if( debug ) std::cout << "BoostFill INFO! H0 tau21 = " << tau21_fj1 << " H1 tau21 = "<< tau21_fj2 << std::endl;
  if( debug ) std::cout << "BoostFill INFO! H0 C2 = " << C2_fj1 << " H1 C2 = "<< C2_fj2 << std::endl;
  if( debug ) std::cout << "BoostFill INFO! H0 D2 = " << D2_fj1 << " H1 D2 = "<< D2_fj2 << std::endl;

  FillHistogram("split12_fj1" + suffix, weight, split12_fj1);
  FillHistogram("split12_fj2" + suffix, weight, split12_fj2);

  FillHistogram("tau21_fj1" + suffix, weight, tau21_fj1);
  FillHistogram("tau21_fj2" + suffix, weight, tau21_fj2);

  FillHistogram("C2_fj1" + suffix, weight, C2_fj1);
  FillHistogram("C2_fj2" + suffix, weight, C2_fj2);

  FillHistogram("D2_fj1" + suffix, weight, D2_fj1);
  FillHistogram("D2_fj2" + suffix, weight, D2_fj2);
  
  if( debug ) std::cout << "BoostFill INFO: Finished BoostFill" << std::endl;
}
