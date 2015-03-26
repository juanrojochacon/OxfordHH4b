// oxford_res_fr.cc

#include "oxford_comb.h"
#include "utils.h"
#include "settings.h"

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

OxfordCombAnalysis::OxfordCombAnalysis(std::string const& sampleName):
Analysis("oxford_comb", sampleName)
{
// ********************* Histogram settings******************

  const double DeltaRmin = 0;
  const double DeltaRmax = 5;

  const double DeltaPhimin = -3;
  const double DeltaPhimax = 3;

  const double DeltaEtamin = -2.5;
  const double DeltaEtamax = 2.5;


  // 2D histograms
  const size_t nbins = 30; const size_t ptmin = 0;  const size_t ptmax = 900;
  const size_t m_min = 0;  const size_t m_max = 180;
 
  // Resolved histograms ************************************************

    // Higgs histograms
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptDijet");
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptDijet1");
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptDijet2");
  BookHistogram(new YODA::Histo1D(30, 0, 500), "ptDijetDijet");

  // Histograms of dijet systems
  BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_Dijetbb");
  BookHistogram(new YODA::Histo1D(20, DeltaPhimin, DeltaPhimax), "DeltaPhi_Dijetbb");
  BookHistogram(new YODA::Histo1D(20, DeltaEtamin, DeltaEtamax), "DeltaEta_Dijetbb");

  BookHistogram(new YODA::Histo1D(20, 0, 200), "mDijet");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "mDijet1");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "mDijet2");

  // 4b system histograms
  BookHistogram(new YODA::Histo1D(20, 200, 1200), "mDijetDijet");
  BookHistogram(new YODA::Histo1D(20, -2.5, 2.5), "yDijetDijet");

  // Dijet distance
  BookHistogram(new YODA::Histo1D(20, -3, 3), "DijetDijet_deltaPhi");

  // 2-D histograms
  BookHistogram(new YODA::Histo2D(nbins, ptmin, ptmax, nbins, ptmin, ptmax), "ptDijetptDijet");
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mDijetmDijet");

  // *********Boosted Histograms ********************************************

  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptFatJet");
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptFatJet1");
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptFatJet2");
  BookHistogram(new YODA::Histo1D(30, 0, 500), "ptFatJetFatJet");

    // 2-D pt histogram
  BookHistogram(new YODA::Histo2D(nbins, ptmin, ptmax, nbins, ptmin, ptmax), "ptFatjetptFatjet");
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mFatjetmFatjet");

  // ********************************************************************


  const std::string tupleSpec = "# signal source weight";
  outputNTuple<<tupleSpec<<std::endl;
}

void OxfordCombAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{

  Analysis::Analyse(signal, weightnorm, fs);

  // Analyse!
  AnalyseResolved(signal, weightnorm, fs);
  AnalyseBoosted(signal, weightnorm, fs);


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
  outputNTuple <<signal <<"\t"<<GetSample()<<"\t"<<weightnorm<<std::endl;
  // Other combinations of kinematical variables could also be useful
  // Need to investigate the kinematics of the 4b final state in more detail
  // Pass event
  Pass(weightnorm);

}


void OxfordCombAnalysis::AnalyseResolved(bool const& signal, double const& weight_norm, finalState const& fs)
{

  // Resolved weight
  double res_weight = weight_norm;

  static double const ResJetR=0.6;
  fastjet::JetDefinition akt(fastjet::antikt_algorithm, ResJetR);

  // Get all the jets with a pT > 40 GeV
  fastjet::ClusterSequence cs_akt(fs, akt);
  std::vector<fastjet::PseudoJet> pt40Jets = sorted_by_pt( cs_akt.inclusive_jets( 0 )  );

  // Not enough jets
  if (pt40Jets.size() < 4 ) return;

  // b-tagging: Loop over the 4 hardest jets in event only (100% efficiency here)
  std::vector<fastjet::PseudoJet> bjets;
  for(int ijet=0; ijet<pt40Jets.size();ijet++) 
  {
    int bQuarks = BTagging(pt40Jets[ijet]);

    if( bQuarks > 0 )
      bjets.push_back(pt40Jets.at(ijet));      
  }

  // We require 4 b-jets in the event, else discard event
  if((int)bjets.size() < 4) 
    return;

  // The next step is to apply the dijet selection
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
  std::vector<fastjet::PseudoJet> res_higgs;
  res_higgs.push_back(bjets.at( jet1_id1) + bjets.at( jet1_id2)); 
  res_higgs.push_back(bjets.at( jet2_id1) + bjets.at( jet2_id2));
  const fastjet::PseudoJet res_dihiggs= res_higgs[0]+res_higgs[1];
  // Sort by higgs pt
  res_higgs = sorted_by_pt(res_higgs);


  // 40 GeV higgs window
  if (abs(res_higgs[0].m() -125) > 40)
    return;
  if (abs(res_higgs[1].m() -125) > 40)
    return;


  // Eta/pT cuts
  const double resEta = 2.5;
  const double respT = 40;
  if (abs(res_higgs[0].eta()) > resEta || res_higgs[0].pt() < respT)
    return;
  if (abs(res_higgs[1].eta()) > resEta || res_higgs[1].pt() < respT)
    return;

  FillHistogram("ptDijet", res_weight, res_higgs[0].pt() );
  FillHistogram("ptDijet", res_weight, res_higgs[1].pt() );

  FillHistogram("ptDijet1", res_weight, res_higgs[0].pt() );
  FillHistogram("ptDijet2", res_weight, res_higgs[1].pt() );

  FillHistogram("ptDijetDijet", res_weight, res_higgs[0].pt() + res_higgs[1].pt());

  FillHistogram("DeltaR_Dijetbb", res_weight, bjets.at(jet1_id1).delta_R(bjets.at(jet1_id2)) );
  FillHistogram("DeltaR_Dijetbb", res_weight, bjets.at(jet2_id1).delta_R(bjets.at(jet2_id2)) );

  FillHistogram("DeltaPhi_Dijetbb", res_weight, bjets.at(jet1_id1).delta_phi_to(bjets.at(jet1_id2)) );
  FillHistogram("DeltaPhi_Dijetbb", res_weight, bjets.at(jet1_id1).delta_phi_to(bjets.at(jet1_id2)) );

  FillHistogram("DeltaEta_Dijetbb", res_weight, bjets.at(jet1_id1).eta() - bjets.at(jet1_id2).eta() );
  FillHistogram("DeltaEta_Dijetbb", res_weight, bjets.at(jet2_id1).eta() - bjets.at(jet2_id2).eta() );

  FillHistogram("mDijet", res_weight, res_higgs[0].m() );
  FillHistogram("mDijet", res_weight, res_higgs[1].m() );

  FillHistogram("mDijet1", res_weight, res_higgs[0].m() );
  FillHistogram("mDijet2", res_weight, res_higgs[1].m() );

  FillHistogram("mDijetDijet", res_weight, res_dihiggs.m() );
  FillHistogram("yDijetDijet", res_weight, res_dihiggs.rapidity() );

    // 2-D Histogram
  FillHistogram("ptDijetptDijet", event_weight, res_higgs[0].pt(), res_higgs[1].pt());
  FillHistogram("mDijetmDijet", event_weight, res_higgs[0].m(), res_higgs[1].m());

}

void OxfordCombAnalysis::AnalyseBoosted(bool const& signal, double const& weight_norm, finalState const& fs)
{
  double boost_weight = weight_norm;
  // Perform jet clustering with anti-kT
  double const boostR=1.0; // To avoid overlapping b's as much as possible
  fastjet::JetDefinition boost_akt(fastjet::antikt_algorithm, boostR);

  // Get all the jets with a pT > 200
  fastjet::ClusterSequence boost_cs_akt(fs, boost_akt);
  std::vector<fastjet::PseudoJet> fatjets = sorted_by_pt( boost_cs_akt.inclusive_jets( 200.00 ) );

  // Not enough jets
  if (fatjets.size() < 2 )
    return;

  // b-tagging: Loop over the 2 hardest jets in event only (100% efficiency here)
  std::vector<fastjet::PseudoJet> fatbjets_unsort;
  for(int ijet=0; ijet<2;ijet++) 
  {
    int bQuarks = BTagging(fatjets[ijet]);

    if( bQuarks > 0 )
      fatbjets_unsort.push_back(fatjets[ijet]);      
  }

  // We require 2 fatjets
  std::vector<fastjet::PseudoJet> fatbjets = sorted_by_pt( fatbjets_unsort );
  if((int)fatbjets.size() == 2)
  {
    // Check if it's a truth higgs
    if (abs(fatbjets[0].m() -125) > 40.0)
      boost_weight = 0;
    if (abs(fatbjets[1].m() -125) > 40.0)
      boost_weight = 0;

    const fastjet::PseudoJet dihiggs= fatbjets[0]+fatbjets[1];
    FillHistogram("ptFatJet", boost_weight, fatbjets[0].pt() );
    FillHistogram("ptFatJet", boost_weight, fatbjets[1].pt() );
    FillHistogram("ptFatJet1", boost_weight, fatbjets[0].pt() );
    FillHistogram("ptFatJet2", boost_weight, fatbjets[1].pt() );
    FillHistogram("ptFatJetFatJet", boost_weight, fatbjets[0].pt() + fatbjets[1].pt());

    // 2-D Histogram
    FillHistogram("ptFatjetptFatjet", event_weight, res_higgs[0].pt(), res_higgs[1].pt());
    FillHistogram("mFatjetmFatjet", event_weight, fatbjets[0].m(), fatbjets[1].m());
  }

}


// ----------------------------------------------------------------------------------

int OxfordCombAnalysis::BTagging( fastjet::PseudoJet const& jet ) const
{
  // Cuts for the b-quark candidates for b-tagging
  double const pt_btagging=5.;

  // Get the jet constituents
  const std::vector<fastjet::PseudoJet>& jet_constituents = jet.constituents();

  int countB=0;
  
  // Loop over constituents and look for b quarks
  // also b quarks must be above some minimum pt
  for(size_t i=0; i<jet_constituents.size(); i++)
  {
    // Flavour of jet constituent
    const int userid= jet_constituents.at(i).user_index();
    const double pt_bcandidate = jet_constituents.at(i).pt();

    if(abs(userid) ==5 )
      if( pt_bcandidate > pt_btagging)
          countB++;
  }
  
  return countB;
}


int OxfordCombAnalysis::CTagging( fastjet::PseudoJet const& jet ) const
{
  // Cuts for the c-quark candidates for c-(mis)tagging
  double const pt_ctagging=5.;

  // Get the jet constituents
  const std::vector<fastjet::PseudoJet>& jet_constituents = jet.constituents();

  int countC=0;
  
  // Loop over constituents and look for c quarks
  // also c quarks must be above some minimum pt
  for(size_t i=0; i<jet_constituents.size(); i++)
  {
    // Flavour of jet constituent
    const int userid= jet_constituents.at(i).user_index();
    const double pt_ccandidate = jet_constituents.at(i).pt();

    if(abs(userid) ==4 )
      if( pt_ccandidate > pt_ctagging)
          countC++;
  }
  
  return countC;
}
