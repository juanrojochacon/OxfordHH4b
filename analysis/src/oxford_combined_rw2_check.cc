// oxford_res_fr.cc

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

#include "oxford_combined_rw2_check.h"
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

// Analysis settings
const int nAnalysis = 3;  const int nCuts = 7;
const std::string aString[nAnalysis] = {"_res", "_inter", "_boost"};
const std::string cString[nCuts] = {"_C0", "_C1a", "_C1b", "_C1c", "_C1d", "_C1e", "_C2"};


OxfordCombinedCheckAnalysis::OxfordCombinedCheckAnalysis(std::string const& sampleName):
Analysis("oxford_combined_check", sampleName)
{
  
  // ********************* Histogram definitions ******************

  for (int i=0; i< nAnalysis; i++)
  {
    BookHistogram(new YODA::Histo1D( nCuts, 0, nCuts ), "CF" + aString[i]);
    BookHistogram(new YODA::Histo1D( nCuts, 0, nCuts ), "CFN" + aString[i]);
  }

}

void OxfordCombinedCheckAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
  Analysis::Analyse(signal, weightnorm, fs);

  bool selected = false;

  // Set initial weight
  const double event_weight = weightnorm;

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

    if (largeRJets_pTcut.size() >= 2) // pT cut
    {
      FillHistogram("CF_boost", event_weight, 2.1);
      FillHistogram("CFN_boost", 1., 2.1);

      if (largeRJets_etacut.size() >= 2) // Eta cut
      {
        FillHistogram("CF_boost", event_weight, 3.1);
        FillHistogram("CFN_boost", 1., 3.1);

        if( largeRJets.size() >= 2 ) // MDT
        {
          FillHistogram("CF_boost", event_weight, 4.1);
          FillHistogram("CFN_boost", 1., 4.1);

          // Higgs mass-window
          const double diffHiggs_0 = fabs(largeRJets[0].m() - 125.);
          const double diffHiggs_1 = fabs(largeRJets[1].m() - 125.);

          if( (diffHiggs_0 < massWindow) && (diffHiggs_1 < massWindow) )
          {
            FillHistogram("CF_boost", event_weight, 5.1);
            FillHistogram("CFN_boost", 1., 5.1);

            // b-tagging weights
            const int nB = nBSubjetsLR_vec[0] + nBSubjetsLR_vec[1];      // Number of true b-subjets
            const int nC = nCSubjetsLR_vec[0] + nCSubjetsLR_vec[1];      // Number of fake b-subjets
            const int nL = nLSubjetsLR_vec[0] + nLSubjetsLR_vec[1];      // Number of fake b-subjets

            // Selection probability
            if (nB+nC+nL == 4) // Selected 4 candidates
            {
              const double P_select_boost = btagProb(4,nB,nC,nL);
              const double dice = ((double) rand() / (double)(RAND_MAX));

              if (dice < P_select_boost && !selected)
              {
                FillHistogram("CF_boost", event_weight, 6.1);
                FillHistogram("CFN_boost", 1., 6.1);
                Pass(event_weight); selected=true;
              }
            }
          }
        }
      }
    }
  }

  // ************************************* Intermediate analysis ********************************************

  if (smallRJets_noCut.size() >= 2 &&  largeRJets_noCut.size() == 1 ) // Clustering
  {
    FillHistogram("CF_inter", event_weight, 1.1);
    FillHistogram("CFN_inter", 1., 1.1);
  }

  if (smallRJets_pTcut.size() >= 2 && largeRJets_pTcut.size() == 1) // pT cut
  {
    FillHistogram("CF_inter", event_weight, 2.1);
    FillHistogram("CFN_inter", 1., 2.1);
  }

  if (smallRJets_etacut.size() >= 2 && largeRJets_etacut.size() == 1) // Eta cut
  {
    FillHistogram("CF_inter", event_weight, 3.1);
    FillHistogram("CFN_inter", 1., 3.1);
  }

  if( smallRJets.size() >= 2 &&  largeRJets.size() == 1 ) // MDT + reco cut
  {
    // Reconstruct Higgs candidates from large-R and small-R jets
    std::vector<fastjet::PseudoJet> higgs_inter;
    std::vector<btagType>  btag_selected_vec;
    fastjet::PseudoJet res_leading_subjet; fastjet::PseudoJet res_subleading_subjet; 

    const bool isRecoInter = Reco_Intermediate( smallRJets, largeRJets[0],
                                                tagType_SR, btag_selected_vec,
                                                res_leading_subjet, res_subleading_subjet,
                                                higgs_inter );

    if( isRecoInter )
    {   
      FillHistogram("CF_inter", event_weight, 4.1);
      FillHistogram("CFN_inter", 1., 4.1);

      // Higgs mass-window
      const double diffHiggs_0 = fabs(higgs_inter[0].m() - 125.);
      const double diffHiggs_1 = fabs(higgs_inter[1].m() - 125.);

      if( (diffHiggs_0 < massWindow) && (diffHiggs_1 < massWindow) )
      {

        FillHistogram("CF_inter", event_weight, 5.1);
        FillHistogram("CFN_inter", 1., 5.1);

        // b-tagging weights
        const int nB = nBSubjetsLR_vec[0] + (btag_selected_vec[0] == BTAG) + (btag_selected_vec[1] == BTAG);      // Number of true b-subjets
        const int nC = nCSubjetsLR_vec[0] + (btag_selected_vec[0] == CTAG) + (btag_selected_vec[1] == CTAG);      // Number of fake b-subjets
        const int nL = nLSubjetsLR_vec[0] + (btag_selected_vec[0] == LTAG) + (btag_selected_vec[1] == LTAG);      // Number of fake b-subjets

        if (nB+nC+nL == 4)
        {
          // Selection probability
          const double P_select_inter = btagProb(4,nB,nC,nL);
          const double dice = ((double) rand() / (double)(RAND_MAX));

          if (dice < P_select_inter && !selected)
          {
            FillHistogram("CF_inter", event_weight, 6.1);
            FillHistogram("CFN_inter", 1., 6.1);
            Pass(event_weight);
            selected =true;
          }
        }
      }
    }
  }


  // ************************************* Resolved analysis ********************************************

  if (smallRJets_noCut.size() >= 4) // Clustering
  {
    FillHistogram("CF_res", event_weight, 1.1);
    FillHistogram("CFN_res", 1., 1.1);

    if (smallRJets_pTcut.size() >= 4) // pT cut
    {
      FillHistogram("CF_res", event_weight, 2.1);
      FillHistogram("CFN_res", 1., 2.1);

      if (smallRJets_etacut.size() >= 4) // Eta cut
      {
        FillHistogram("CF_res", event_weight, 3.1);
        FillHistogram("CFN_res", 1., 3.1);

        if( smallRJets.size() >= 4 )
        {
          FillHistogram("CF_res", event_weight, 4.1);
          FillHistogram("CFN_res", 1., 4.1);

          // Reconstruct Higgs candidates from small-R jets
          std::vector<fastjet::PseudoJet> higgs_res;
          std::vector<fastjet::PseudoJet> higgs0_res;
          std::vector<fastjet::PseudoJet> higgs1_res;

          Reco_Resolved( smallRJets, higgs_res, higgs0_res, higgs1_res );

          // Higgs mass-window
          const double diffHiggs_0 = fabs(higgs_res[0].m() - 125.);
          const double diffHiggs_1 = fabs(higgs_res[1].m() - 125.);

          if( ( diffHiggs_0 < massWindow ) && ( diffHiggs_1 < massWindow ) )
          {
            FillHistogram("CF_res", event_weight, 5.1);
            FillHistogram("CFN_res", 1., 5.1);

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
              const double P_select_resol = btagProb( 4, nB, nC, nL);
              const double dice = ((double) rand() / (double)(RAND_MAX));

              if (dice < P_select_resol && !selected)
              {
                FillHistogram("CF_res", event_weight, 6.1);
                FillHistogram("CFN_res", 1., 6.1);
                Pass(event_weight);
                selected =true;
              }
            }
          }
        }
      }
    }
  }

  if (!selected)
    Cut ("Uncategorised", event_weight);

  return;
}

// Small-R B-tagging
void OxfordCombinedCheckAnalysis::BTagging( std::vector<fastjet::PseudoJet> const& jets_vec, std::vector<btagType>& btag_vec  )
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

void OxfordCombinedCheckAnalysis::BTagging( std::vector<fastjet::PseudoJet> const& largeRJets,
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


void OxfordCombinedCheckAnalysis::Reco_Resolved( std::vector<fastjet::PseudoJet> const& bjets, // Input b-jets
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


bool OxfordCombinedCheckAnalysis::Reco_Intermediate(  std::vector<fastjet::PseudoJet> const& bjets, 
                                                    fastjet::PseudoJet const& fatjet, 
                                                    std::vector<btagType> const& btag_vec,
                                                    std::vector<btagType>&  btag_selected_vec,
                                                    fastjet::PseudoJet& lead_subjet,
                                                    fastjet::PseudoJet& sublead_subjet,
                                                    std::vector<fastjet::PseudoJet>& higgs_vec )
{
  // Identify small-R jets separated from merged Higgs
  std::vector<fastjet::PseudoJet> bjets_separated;
  for( size_t i = 0; i < bjets.size(); i++ )
  {  
    double dR = bjets.at(i).delta_R(fatjet);
    if( dR < 1.2 ) continue;
    
    if (bjets_separated.size() < 2)
    {
      bjets_separated.push_back( bjets[i] );
      btag_selected_vec.push_back(btag_vec[i]);
    }
  }
  
  if( bjets_separated.size() < 2 ) return false;
  
  // Get the pairing that minimizes |m_dj1 - m_dj2|
  double mdj_diff_min = std::numeric_limits<double>::infinity();
  int bjet_id1(-1), bjet_id2(-1);
  for(int ijet=0; ijet < 2; ijet++)
    for(int jjet=0; jjet < 2; jjet++)
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
