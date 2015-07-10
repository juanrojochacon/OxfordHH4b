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
static const bool debug = false;


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
  
  const double m_HH_min = 0.;
  const double m_HH_max = 600.; 
  
  const double pt_HH_min = 0.;
  const double pt_HH_max = 300.;
  
  const int nbins = 30;
  
  // ********************* Histogram definitions ******************

  const int nAnalysis = 3;  const int nCuts = 8;
  const std::string aString[nAnalysis] = {"_res", "_inter", "_boost"};
  const std::string cString[nCuts] = {"_C0", "_C1a", "_C1b", "_C1c", "_C1d", "_C1e", "_C2", "_C3"};

  for (int i=0; i< nAnalysis; i++)
  {
    BookHistogram(new YODA::Histo1D( nCuts+1, 0, nCuts+1 ), "CF" + aString[i]);
    BookHistogram(new YODA::Histo1D( nCuts+1, 0, nCuts+1 ), "CFN" + aString[i]);

    for (int j=0; j< nCuts; j++)
    {
      const std::string suffix = aString[i] + cString[j];

      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1" + suffix);

      BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1" + suffix);

      BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH" + suffix);

      BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH" + suffix);
      BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH" + suffix);

      BookHistogram(new YODA::Histo2D(nbins, pt_min, pt_max, nbins, pt_min, pt_max), "ptHptH" + suffix);
      BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mHmH" + suffix);

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

  const std::string tupleSpec = "# signal source weight pt_H0 pt_H1 pt_HH m_H0 m_H1 m_HH dR_HH dPhi_HH dEta_HH";
  outputNTuple<<tupleSpec<<std::endl;

  const std::string root = "." + GetRoot() + GetSample() + "/";

  const std::string resDir = root+"resNTuple.dat";
  const std::string intDir = root+"intNTuple.dat";
  const std::string bstDir = root+"bstNTuple.dat";

  resNTuple.open(resDir.c_str());
  intNTuple.open(intDir.c_str());
  bstNTuple.open(bstDir.c_str());

  resNTuple << tupleSpec <<" pt_H0_sub0 pt_H0_sub1 pt_H1_sub0 pt_H1_sub1"<<std::endl;
  intNTuple << tupleSpec <<" split12_fj tau21_fj C2_fj D2_fj"<<std::endl;
  bstNTuple << tupleSpec <<" split12_fj1 split12_fj2 tau21_fj1 tau21_fj2 C2_fj1 C2_fj2 D2_fj1 D2_fj2"<<std::endl;

  // ********************* Cutflow histograms  **********************

  // Category overlap
  BookHistogram(new YODA::Histo1D( 8, 0, 8 ), "Categories_C1a");
  BookHistogram(new YODA::Histo1D( 8, 0, 8 ), "Categories_C1b");
  BookHistogram(new YODA::Histo1D( 8, 0, 8 ), "Categories_C1c");
  BookHistogram(new YODA::Histo1D( 8, 0, 8 ), "Categories_C1d");
  BookHistogram(new YODA::Histo1D( 8, 0, 8 ), "Categories_C1e");
  BookHistogram(new YODA::Histo1D( 8, 0, 8 ), "Categories_C2");
  BookHistogram(new YODA::Histo1D( 8, 0, 8 ), "Categories_C3");

}

void OxfordCombinedRW2Analysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
  Analysis::Analyse(signal, weightnorm, fs);

  // Set initial weight
  const double event_weight = weightnorm;
  bool selected = false;

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
    if (smallRJets_pTcut[i].eta() <= SR_maxEta)
      smallRJets_etacut.push_back(smallRJets_pTcut[i]);

  // Cluster large-R jets
  const fastjet::JetDefinition akt_boost(fastjet::antikt_algorithm, BoostJetR);
  const fastjet::ClusterSequence cs_akt_bst(fs, akt_boost);
  const std::vector<fastjet::PseudoJet> largeRJets_noCut = sorted_by_pt( cs_akt_bst.inclusive_jets()  ); 
  const std::vector<fastjet::PseudoJet> largeRJets_pTcut = sorted_by_pt( cs_akt_bst.inclusive_jets( LR_minPT ) ); 

  // Eta cut
  std::vector<fastjet::PseudoJet> largeRJets_etacut;
  for (size_t i=0; i<largeRJets_pTcut.size(); i++)
    if (largeRJets_pTcut[i].eta() <= LR_maxEta)
      largeRJets_etacut.push_back(largeRJets_pTcut[i]);

  // Cluster small-R track jets 
  const fastjet::JetDefinition jd_subjets(fastjet::antikt_algorithm, GAjetR);
  const fastjet::ClusterSequence cs_subjets(fs, jd_subjets);
  const std::vector<fastjet::PseudoJet> trackjets_nocut = sorted_by_pt( cs_subjets.inclusive_jets()  );
  const std::vector<fastjet::PseudoJet> trackjets_pTcut = sorted_by_pt( cs_subjets.inclusive_jets( TJ_minPT ) );

  // Eta cut
  std::vector<fastjet::PseudoJet> trackJets_etacut;
  for (size_t i=0; i<trackjets_pTcut.size(); i++)
    if (trackjets_pTcut[i].eta() <= TJ_maxEta)
      trackJets_etacut.push_back(trackjets_pTcut[i]);

  // Final sorted jets
  const std::vector<fastjet::PseudoJet> smallRJets = sorted_by_pt(smallRJets_etacut);
  const std::vector<fastjet::PseudoJet> largeRJets = sorted_by_pt(largeRJets_etacut);
  const std::vector<fastjet::PseudoJet> trackJets = sorted_by_pt(trackJets_etacut);

  // ***************************************** Initial histograms **********************************************

  FillHistogram("CF_res", event_weight, 0.1);
  FillHistogram("CF_inter", event_weight, 0.1);
  FillHistogram("CF_boost", event_weight, 0.1);

  FillHistogram("CFN_res", 1., 0.1);
  FillHistogram("CFN_inter", 1., 0.1);
  FillHistogram("CFN_boost", 1., 0.1);

  // Resolved initial histograms
  if (smallRJets_noCut.size() >= 4)
  {
    FillHistogram("CF_res", event_weight, 1.1);
    FillHistogram("CFN_res", 1., 1.1);

    if (smallRJets_pTcut.size() >= 4)
    {
      FillHistogram("CF_res", event_weight, 2.1);
      FillHistogram("CFN_res", 1., 2.1);
    }
  }

  // Boosted initial histograms
  if (largeRJets_noCut.size() >= 2)
  {
    FillHistogram("CF_boost", event_weight, 1.1);
    FillHistogram("CFN_boost", 1., 1.1);

    if (largeRJets_pTcut.size() >= 2)
    {
      FillHistogram("CF_boost", event_weight, 2.1);
      FillHistogram("CFN_boost", 1., 2.1);
    }
  }

  // Intermediate initial histograms
  if (largeRJets_noCut.size() == 1 && smallRJets_noCut.size() >= 2)
  {
    FillHistogram("CF_inter", event_weight, 1.1);
    FillHistogram("CFN_inter", 1., 1.1);

    if (largeRJets_pTcut.size() == 1 && smallRJets_pTcut.size() >= 2)
    {
      FillHistogram("CF_inter", event_weight, 2.1);
      FillHistogram("CFN_inter", 1., 2.1);
    }
  }

  // ***************************************** B-Tagging **********************************************

  // b-tagging for large-R jets
  std::vector<int> nBSubjetsLR_vec;    // Vector specifying how many real b subjets there are
  BTagging( largeRJets, trackJets, nBSubjetsLR_vec);
  if( largeRJets.size() != nBSubjetsLR_vec.size() )
    std::cout << "ERROR: b-tagging vector sizes don't match number of fat jets" << std::endl;

  // b-tagging for small-R jets
  std::vector<bool> isFakeSR_vec;  // Vector specifying if each jet is fake or not
  BTagging( smallRJets, isFakeSR_vec );

  // **************************************** Boosted analysis *********************************************

  if( largeRJets.size() >= 2 ) // Eta cut
  {
    // Third cut flow fill
    HiggsFill(largeRJets[0], largeRJets[1], "boost", 3, event_weight);

    // Check if jets are mass-drop tagged
    const fastjet::JetDefinition CA10(fastjet::cambridge_algorithm, 1.0);
    const fastjet::MassDropTagger md_tagger(mu, ycut);

    // Recluster constituents
    const fastjet::ClusterSequence cs_sub_0( largeRJets[0].constituents(), CA10);
    const fastjet::ClusterSequence cs_sub_1( largeRJets[1].constituents(), CA10);
    const fastjet::PseudoJet ca_jet_0 = sorted_by_pt(cs_sub_0.inclusive_jets())[0];
    const fastjet::PseudoJet ca_jet_1 = sorted_by_pt(cs_sub_1.inclusive_jets())[0];

    const fastjet::PseudoJet tagged_jet_0 = md_tagger(ca_jet_0);
    const fastjet::PseudoJet tagged_jet_1 = md_tagger(ca_jet_1);

    // Mass-drop tagged
    if ( tagged_jet_0 != 0 && tagged_jet_1 != 0 )
    {
      HiggsFill(largeRJets[0], largeRJets[1], "boost", 4, event_weight);
      BoostFill(largeRJets[0], largeRJets[1], "boost", 4, event_weight);

      // Higgs mass-window
      const double diffHiggs_0 = fabs(largeRJets[0].m() - 125.);
      const double diffHiggs_1 = fabs(largeRJets[1].m() - 125.);

      if( (diffHiggs_0 < massWindow) && (diffHiggs_1 < massWindow) )
      {
        HiggsFill(largeRJets[0], largeRJets[1], "boost", 5, event_weight);
        BoostFill(largeRJets[0], largeRJets[1], "boost", 5, event_weight);

        // b-tagging weights
        const double nB = nBSubjetsLR_vec[0] + nBSubjetsLR_vec[1];      // Number of true b-subjets
        const double nF = 4 - nB;  // Number of fake b-subjets

        // Reweighted event weight
        const double boost_weight = pow(btag_prob,nB)*pow(btag_mistag,nF)*event_weight;
        const fastjet::PseudoJet dihiggs_boost = largeRJets[0] + largeRJets[1];

        HiggsFill(largeRJets[0], largeRJets[1], "boost", 6, boost_weight);
        BoostFill(largeRJets[0], largeRJets[1], "boost", 6, boost_weight);

        // Final exclusive booking
        if (!selected)
        {
          selected = true;

          HiggsFill(largeRJets[0], largeRJets[1], "boost", 7, boost_weight);
          BoostFill(largeRJets[0], largeRJets[1], "boost", 7, boost_weight);

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

  // ************************************* Intermediate analysis ********************************************

  // Intermediate
  if( smallRJets.size() >= 2 &&  largeRJets.size() == 1 )
  {
    // Reconstruct Higgs candidates from large-R and small-R jets
    std::vector<fastjet::PseudoJet> higgs_inter; int nBJets_SR = 0;
    const bool isRecoInter = Reco_Intermediate( smallRJets, isFakeSR_vec, largeRJets[0], nBJets_SR, higgs_inter );

    if( isRecoInter )
    {
      HiggsFill(higgs_inter[0], higgs_inter[1], "inter", 3, event_weight);

      // Check if large-R jet is mass-drop tagged
      const fastjet::JetDefinition CA10(fastjet::cambridge_algorithm, 1.0);
      const fastjet::MassDropTagger md_tagger(mu, ycut);

      // Recluster constituents
      const fastjet::ClusterSequence cs_sub( largeRJets[0].constituents(), CA10);
      const fastjet::PseudoJet ca_jet = sorted_by_pt(cs_sub.inclusive_jets())[0];
      const fastjet::PseudoJet tagged_jet = md_tagger(ca_jet);

      // Mass-drop tagged
      if ( tagged_jet != 0 )
      {
        HiggsFill(higgs_inter[0], higgs_inter[1], "inter", 4, event_weight);
        BoostFill(largeRJets[0], "inter", 4, event_weight);

        // Higgs mass-window
        const double diffHiggs_0 = fabs(higgs_inter[0].m() - 125.);
        const double diffHiggs_1 = fabs(higgs_inter[1].m() - 125.);

        if( (diffHiggs_0 < massWindow) && (diffHiggs_1 < massWindow) )
        {
          HiggsFill(higgs_inter[0], higgs_inter[1], "inter", 5, event_weight);
          BoostFill(largeRJets[0], "inter", 5, event_weight);

          // Determine number of fake bJets
          const int nB = nBJets_SR + nBSubjetsLR_vec[0];
          const int nF = 4 - nB;

          // Error checking
          if (nB > 4)
          {
            std::cerr << "ERROR: number of reconstructed b quarks > 4!"<<std::endl;
            exit(-1);
          }

          // Reweighted event weight
          const double inter_weight = pow(btag_prob,nB)*pow(btag_mistag,nF)*event_weight;
          const fastjet::PseudoJet dihiggs_inter = higgs_inter[0] + higgs_inter[1];

          HiggsFill(higgs_inter[0], higgs_inter[1], "inter", 6, inter_weight);
          BoostFill(largeRJets[0], "inter", 6, inter_weight);

          if (!selected)
          {
            selected = true;

            // Exclusivity cut
            HiggsFill(higgs_inter[0], higgs_inter[1], "inter", 7, inter_weight);
            BoostFill(largeRJets[0], "inter", 7, inter_weight);

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
                      << split12 << "\t"
                      << tau21 << "\t"
                      << C2 << "\t"
                      << D2 << "\t"
                      <<std::endl;


            // Final
            Pass(inter_weight);
            Cut("IntermediateCut", event_weight - inter_weight );
          }
        }
      }
    }
  }

  // ************************************* Resolved analysis ********************************************

  // Resolved
  if( smallRJets.size() >= 4 )
  {
    // Reconstruct Higgs candidates from small-R jets
    std::vector<fastjet::PseudoJet> higgs_res;
    std::vector<fastjet::PseudoJet> higgs0_res;
    std::vector<fastjet::PseudoJet> higgs1_res;

    Reco_Resolved( smallRJets, higgs_res, higgs0_res, higgs1_res );
   
    HiggsFill( higgs_res[0], higgs_res[1], "res", 3, event_weight );
    HiggsFill( higgs_res[0], higgs_res[1], "res", 4, event_weight ); // no MDT

    // Higgs mass-window
    const double diffHiggs_0 = fabs(higgs_res[0].m() - 125.);
    const double diffHiggs_1 = fabs(higgs_res[1].m() - 125.);

    if( (diffHiggs_0 < massWindow) && (diffHiggs_1 < massWindow) )
    {
      HiggsFill( higgs_res[0], higgs_res[1], "res", 5, event_weight );

      // Determine number of real and fake b-jets
      int nB = 0;
      for (int i=0; i<4; i++)
        if (!isFakeSR_vec[i])
          nB++;

      const int nF = 4 - nB;

      // Reweighted event weight
      const double res_weight = pow(btag_prob,nB)*pow(btag_mistag,nF)*event_weight;
      const fastjet::PseudoJet dihiggs_res = higgs_res[0] + higgs_res[1];

      HiggsFill(higgs_res[0], higgs_res[1], "res", 6, res_weight);

      if (!selected)
      {
        selected = true;

        HiggsFill(higgs_res[0], higgs_res[1], "res", 7, res_weight);
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

  // ************************************* Remaining ********************************************
 
  /*
  
  if( isRes_C1a && !isInter_C1a && !isBoost_C1a )	FillHistogram("Categories_C1a", 1.0, 0.1);
  else if( !isRes_C1a && isInter_C1a && !isBoost_C1a )	FillHistogram("Categories_C1a", 1.0, 1.1);
  else if( !isRes_C1a && !isInter_C1a && isBoost_C1a )	FillHistogram("Categories_C1a", 1.0, 2.1);
  else if( isRes_C1a && isInter_C1a && !isBoost_C1a )	FillHistogram("Categories_C1a", 1.0, 3.1);
  else if( isRes_C1a && !isInter_C1a && isBoost_C1a )	FillHistogram("Categories_C1a", 1.0, 4.1);
  else if( !isRes_C1a && isInter_C1a && isBoost_C1a )	FillHistogram("Categories_C1a", 1.0, 5.1);
  else if( isRes_C1a && isInter_C1a && isBoost_C1a )	FillHistogram("Categories_C1a", 1.0, 6.1);

  //===================================================
  // C1b: Jet pt cuts
  //===================================================

   
  if( isRes_C1b && !isInter_C1b && !isBoost_C1b )	FillHistogram("Categories_C1b", 1.0, 0.1);
  else if( !isRes_C1b && isInter_C1b && !isBoost_C1b )	FillHistogram("Categories_C1b", 1.0, 1.1);
  else if( !isRes_C1b && !isInter_C1b && isBoost_C1b )	FillHistogram("Categories_C1b", 1.0, 2.1);
  else if( isRes_C1b && isInter_C1b && !isBoost_C1b )	FillHistogram("Categories_C1b", 1.0, 3.1);
  else if( isRes_C1b && !isInter_C1b && isBoost_C1b )	FillHistogram("Categories_C1b", 1.0, 4.1);
  else if( !isRes_C1b && isInter_C1b && isBoost_C1b )	FillHistogram("Categories_C1b", 1.0, 5.1);
  else if( isRes_C1b && isInter_C1b && isBoost_C1b )	FillHistogram("Categories_C1b", 1.0, 6.1);
 
  //===================================================
  // C1c: Jet pt and eta cuts + Higgs mass window cut
  //===================================================
  
   
  
  if( isRes_C1c && !isInter_C1c && !isBoost_C1c )	FillHistogram("Categories_C1c", 1.0, 0.1);
  else if( !isRes_C1c && isInter_C1c && !isBoost_C1c )	FillHistogram("Categories_C1c", 1.0, 1.1);
  else if( !isRes_C1c && !isInter_C1c && isBoost_C1c )	FillHistogram("Categories_C1c", 1.0, 2.1);
  else if( isRes_C1c && isInter_C1c && !isBoost_C1c )	FillHistogram("Categories_C1c", 1.0, 3.1);
  else if( isRes_C1c && !isInter_C1c && isBoost_C1c )	FillHistogram("Categories_C1c", 1.0, 4.1);
  else if( !isRes_C1c && isInter_C1c && isBoost_C1c )	FillHistogram("Categories_C1c", 1.0, 5.1);
  else if( isRes_C1c && isInter_C1c && isBoost_C1c )	FillHistogram("Categories_C1c", 1.0, 6.1);

  //=========================================
  // C1d: All jet cuts
  //=========================================

  
  if( isRes_C1d && !isInter_C1d && !isBoost_C1d )	FillHistogram("Categories_C1d", 1.0, 0.1);
  else if( !isRes_C1d && isInter_C1d && !isBoost_C1d )	FillHistogram("Categories_C1d", 1.0, 1.1);
  else if( !isRes_C1d && !isInter_C1d && isBoost_C1d )	FillHistogram("Categories_C1d", 1.0, 2.1);
  else if( isRes_C1d && isInter_C1d && !isBoost_C1d )	FillHistogram("Categories_C1d", 1.0, 3.1);
  else if( isRes_C1d && !isInter_C1d && isBoost_C1d )	FillHistogram("Categories_C1d", 1.0, 4.1);
  else if( !isRes_C1d && isInter_C1d && isBoost_C1d )	FillHistogram("Categories_C1d", 1.0, 5.1);
  else if( isRes_C1d && isInter_C1d && isBoost_C1d )	FillHistogram("Categories_C1d", 1.0, 6.1);
  
  //===================================================
  // C1e: All jet cuts + Higgs mass window cut
  //===================================================
  
  


  
  if( isRes_C1e && !isInter_C1e && !isBoost_C1e )	FillHistogram("Categories_C1e", 1.0, 0.1);
  else if( !isRes_C1e && isInter_C1e && !isBoost_C1e )	FillHistogram("Categories_C1e", 1.0, 1.1);
  else if( !isRes_C1e && !isInter_C1e && isBoost_C1e )	FillHistogram("Categories_C1e", 1.0, 2.1);
  else if( isRes_C1e && isInter_C1e && !isBoost_C1e )	FillHistogram("Categories_C1e", 1.0, 3.1);
  else if( isRes_C1e && !isInter_C1e && isBoost_C1e )	FillHistogram("Categories_C1e", 1.0, 4.1);
  else if( !isRes_C1e && isInter_C1e && isBoost_C1e )	FillHistogram("Categories_C1e", 1.0, 5.1);
  else if( isRes_C1e && isInter_C1e && isBoost_C1e )	FillHistogram("Categories_C1e", 1.0, 6.1);
  */
  if (!selected)
    return Cut ("Uncategorised", event_weight);
}


void OxfordCombinedRW2Analysis::BTagging( std::vector<fastjet::PseudoJet> const& jets_vec, std::vector<bool>& isFake_vec  )
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

void OxfordCombinedRW2Analysis::BTagging( std::vector<fastjet::PseudoJet> const& largeRJets, std::vector<fastjet::PseudoJet> const& trackjets, std::vector<int>& nBSubJets_vec )
{
    // Loop over all fat jets
    for( size_t i=0; i<largeRJets.size(); i++)
    {
      // Get ghost associated track jets
      std::vector<fastjet::PseudoJet> subjets;
      get_assoc_trkjets( largeRJets.at(i), trackjets, subjets, false);

      int nBSubJets = 0;

      // Loop over subjets
      const int nSubjets = std::min((int)subjets.size(), 2); // Restrict to leading 2 subjets
      for( int j=0; j < nSubjets; j++)
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
          {
            isBJet = true;
            break;
          }
        }

        if (isBJet)
          nBSubJets++;
           
      }//end of loop over subjets
  	   
      nBSubJets_vec.push_back(nBSubJets);

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
                                                  std::vector<bool> const& isFakeSR_vec,
                                                  fastjet::PseudoJet const& fatjet, 
                                                  int& nBjets,
                                                  std::vector<fastjet::PseudoJet>& higgs_vec )
{
  // Identify small-R jets separated from merged Higgs
  std::vector<fastjet::PseudoJet> bjets_separated;
  std::vector<bool> bjets_separated_isFake;
  
  for( size_t i = 0; i < bjets.size(); i++ )
  {  
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
  
  const std::string cutID[8] = {"0","1a","1b","1c","1d","1e","2","3"};
  const std::string cutStr = "_C" + cutID[cut];
  const std::string suffix = "_" + analysis + cutStr;
    
  // Record cutflow
  FillHistogram("CF_" +analysis, weight, cut + 0.1);
  FillHistogram("CFN_"+analysis, 1., cut + 0.1);

  // Histograms for reconstructed Higgs candidates
  FillHistogram("pt_H0" + suffix, weight, H0.pt());
  FillHistogram("pt_H1" + suffix, weight, H1.pt());

  FillHistogram("m_H0" + suffix, weight, H0.m());
  FillHistogram("m_H1" + suffix, weight, H1.m());

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
  const std::string cutID[8] = {"0","1a","1b","1c","1d","1e","2","3"};
  const std::string cutStr = "_C" + cutID[cut];
  const std::string suffix = "_" + analysis + cutStr;

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

  const std::string cutID[8] = {"0","1a","1b","1c","1d","1e","2","3"};
  const std::string cutStr = "_C" + cutID[cut];
  const std::string suffix = "_" + analysis + cutStr;

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
