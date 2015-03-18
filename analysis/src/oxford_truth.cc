// oxford_res_fr.cc

#include "oxford_truth.h"
#include "utils.h"
#include "settings.h"

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

OxfordTruthAnalysis::OxfordTruthAnalysis(std::string const& sampleName):
Analysis("oxford_truth", sampleName)
{
// ********************* Histogram settings******************

  const double DeltaRmin = 0;
  const double DeltaRmax = 5;

  const double DeltaPhimin = -3;
  const double DeltaPhimax = 3;

  const double DeltaEtamin = -2.5;
  const double DeltaEtamax = 2.5;

  // ********************* Histograms before cuts ***************************

  // Mass cross check
  BookHistogram(new YODA::Histo1D(30, -10, 10), "massCrossCheck");

  // Higgs histograms
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptH");
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptH1");
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptH2");
  BookHistogram(new YODA::Histo1D(30, 0, 500), "ptHH");

  // 2-D pt histogram
  const size_t nbins = 30; const size_t ptmin = 0;  const size_t ptmax = 900;
  BookHistogram(new YODA::Histo2D(nbins, ptmin, ptmax, nbins, ptmin, ptmax), "ptHptH");
  
  // 2-D mass histogram
  const size_t m_min = 0;  const size_t m_max = 200;
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max), "mHmH");
 
  // 2-D deltaR histograms
  BookHistogram(new YODA::Histo2D(nbins, ptmin, ptmax, nbins, DeltaRmin, DeltaRmax), "deltaR_pt_H0");
  BookHistogram(new YODA::Histo2D(nbins, ptmin, ptmax, nbins, DeltaRmin, DeltaRmax), "deltaR_pt_H1");
  BookHistogram(new YODA::Histo2D(nbins, ptmin, ptmax, nbins, DeltaRmin, DeltaRmax), "deltaR_pt_bothH");

  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, DeltaRmin, DeltaRmax), "deltaR_m_H0");
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, DeltaRmin, DeltaRmax), "deltaR_m_H1");
  BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, DeltaRmin, DeltaRmax), "deltaR_m_bothH");
 
  // Histograms of dijet systems
  BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_Hbb");
  BookHistogram(new YODA::Histo1D(20, DeltaPhimin, DeltaPhimax), "DeltaPhi_Hbb");
  BookHistogram(new YODA::Histo1D(20, DeltaEtamin, DeltaEtamax), "DeltaEta_Hbb");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "mH");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "mH1");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "mH2");

  // 4b system histograms
  BookHistogram(new YODA::Histo1D(20, 200, 1200), "mHH");
  BookHistogram(new YODA::Histo1D(20, -2.5, 2.5), "yHH");

  // Dijet distance
  BookHistogram(new YODA::Histo1D(20, -3, 3), "HH_deltaPhi");

  // DeltaR between the Higgs and each b-quark
  BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_H0b0");
  BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_H0b1");
  BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_H1b0");
  BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_H1b1");


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

  // *********Boosted Histograms ********************************************

  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptFatJet");
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptFatJet1");
  BookHistogram(new YODA::Histo1D(30, 0, 900), "ptFatJet2");
  BookHistogram(new YODA::Histo1D(30, 0, 500), "ptFatJetFatJet");

  // ********************************************************************


  const std::string tupleSpec = "# signal source weight";
  outputNTuple<<tupleSpec<<std::endl;
}

void OxfordTruthAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
  // Only for signal for now
  if (!signal) return;
  if (fs.size() != 4) std::cout << "ERROR: Please run truth analysis on unshowered samples!"<<std::endl; 

  Analysis::Analyse(signal, weightnorm, fs);

  // Set initial weight
  const double event_weight = weightnorm;

  std::vector<fastjet::PseudoJet> higgs1bb;
  std::vector<fastjet::PseudoJet> higgs2bb;

  // Get ID of first higgs
  const int h1ID = fs[0].user_info<JetInfo>().motherID;
  for (size_t i=0; i<4; i++)
  {
    if (fs[i].user_info<JetInfo>().motherID == h1ID)
      higgs1bb.push_back(fs[i]);
    else
      higgs2bb.push_back(fs[i]);
  }

  // Sort by pt
  higgs1bb = sorted_by_pt(higgs1bb);
  higgs2bb = sorted_by_pt(higgs2bb);

  std::vector<fastjet::PseudoJet> higgs;
  std::vector<fastjet::PseudoJet> higgs_unsrt;
  higgs_unsrt.push_back(higgs1bb[0] + higgs1bb[1]);
  higgs_unsrt.push_back(higgs2bb[0] + higgs2bb[1]);
  
  // Sort by pt: both Higgses and the corresponding 
  // vectors of b-quarks
  if( higgs_unsrt[0].pt() < higgs_unsrt[1].pt() ){
	higgs.push_back( higgs_unsrt[1] );
	higgs.push_back( higgs_unsrt[0] );
        
        // Swap also the corresponding vectors of b-quarks
        std::vector<fastjet::PseudoJet> placeholder = higgs1bb;
        higgs1bb = higgs2bb;
        higgs2bb = placeholder;
  }
  else{	
	higgs.push_back( higgs_unsrt[0] );
	higgs.push_back( higgs_unsrt[1] );
  }

  const fastjet::PseudoJet diHiggs = higgs[0] + higgs[1];

  // Delta R
  const double deltaR_H0b0 = higgs[0].delta_R(higgs1bb[0]);
  const double deltaR_H0b1 = higgs[0].delta_R(higgs1bb[1]);

  const double deltaR_H1b0 = higgs[1].delta_R(higgs1bb[0]);
  const double deltaR_H1b1 = higgs[1].delta_R(higgs1bb[1]);

  const double deltaR_H0bb = higgs1bb[0].delta_R(higgs1bb[1]);
  const double deltaR_H1bb = higgs2bb[0].delta_R(higgs2bb[1]);
  
  // Mass cross check
  fastjet::PseudoJet b1 = higgs1bb[0];
  fastjet::PseudoJet b2 = higgs1bb[1];
  const double mass_H0_FastJet = higgs[0].m();

  const double mass2_H0_byHand = b1.m2() + b2.m2() + 2 * ( b1.e()*b2.e() - b1.px()*b2.px() - b1.py()*b2.py() - b1.pz()*b2.pz() );
  const double mass_H0_byHand = (mass2_H0_byHand/abs(mass2_H0_byHand)) * sqrt(abs(mass2_H0_byHand));

  std::cout << "FastJet mass " << mass_H0_FastJet << " mass by hand " << mass_H0_byHand <<std::endl;
  std::cout << "b1 mass: " << b1.m() << " b2 mass: " << b2.m() << std::endl;
  // *******************************************************************************
  // Histograms
  
  FillHistogram("massCrossCheck", event_weight, mass_H0_FastJet - mass_H0_byHand );

  FillHistogram("ptH", event_weight, higgs[0].pt() );
  FillHistogram("ptH", event_weight, higgs[1].pt() );

  FillHistogram("ptH1", event_weight, higgs[0].pt() );
  FillHistogram("ptH2", event_weight, higgs[1].pt() );

  FillHistogram("ptHH", event_weight, higgs[0].pt() + higgs[1].pt());

  FillHistogram("DeltaR_Hbb", event_weight, higgs1bb[0].delta_R(higgs1bb[1]) );
  FillHistogram("DeltaR_Hbb", event_weight, higgs2bb[0].delta_R(higgs2bb[1]) );

  FillHistogram("DeltaPhi_Hbb", event_weight, higgs1bb[0].delta_phi_to(higgs1bb[1]) );
  FillHistogram("DeltaPhi_Hbb", event_weight, higgs2bb[0].delta_phi_to(higgs2bb[1]) );

  FillHistogram("DeltaEta_Hbb", event_weight, higgs1bb[0].eta() - higgs1bb[1].eta() );
  FillHistogram("DeltaEta_Hbb", event_weight, higgs2bb[0].eta() - higgs2bb[1].eta() );

  FillHistogram("mH", event_weight, higgs[0].m() );
  FillHistogram("mH", event_weight, higgs[1].m() );
  FillHistogram("mH1", event_weight, higgs[0].m() );
  FillHistogram("mH2", event_weight, higgs[1].m() );

  FillHistogram("mHH", event_weight, diHiggs.m() );
  FillHistogram("yHH", event_weight, diHiggs.rapidity() );

  // DeltaR between Higgs and each decay product
  FillHistogram("DeltaR_H0b0", event_weight, deltaR_H0b0); 
  FillHistogram("DeltaR_H0b1", event_weight, deltaR_H0b1); 
  FillHistogram("DeltaR_H1b0", event_weight, deltaR_H1b0); 
  FillHistogram("DeltaR_H1b1", event_weight, deltaR_H1b1); 

  // 2-D Histogram
  FillHistogram("ptHptH", event_weight, higgs[0].pt(), higgs[1].pt());
  FillHistogram("mHmH", event_weight, higgs[0].m(), higgs[1].m());

  FillHistogram("deltaR_pt_H0", event_weight, higgs[0].pt(), deltaR_H0bb);
  FillHistogram("deltaR_pt_H1", event_weight, higgs[1].pt(), deltaR_H1bb);
  FillHistogram("deltaR_pt_bothH", event_weight, higgs[0].pt(), deltaR_H0bb);
  FillHistogram("deltaR_pt_bothH", event_weight, higgs[1].pt(), deltaR_H1bb);

  FillHistogram("deltaR_m_H0", event_weight, higgs[0].m(), deltaR_H0bb);
  FillHistogram("deltaR_m_H1", event_weight, higgs[1].m(), deltaR_H1bb);
  FillHistogram("deltaR_m_bothH", event_weight, higgs[0].m(), deltaR_H0bb);
  FillHistogram("deltaR_m_bothH", event_weight, higgs[1].m(), deltaR_H1bb);
  
  // ********************************* Resolved efficiencies ***************************************

  // Resolved weight
  double res_weight = event_weight;

  static double const ResJetR=0.4; // To avoid overlapping b's as much as possible
  fastjet::JetDefinition akt(fastjet::antikt_algorithm, ResJetR);

  // Cluster all particles
  // The cluster sequence has to be saved to be used for jet substructure
  fastjet::ClusterSequence cs_akt(fs, akt);
  // Get all the jets (no pt cut here)
  std::vector<fastjet::PseudoJet> bjets = sorted_by_pt( cs_akt.inclusive_jets()  );

  // We require at least 4 jets in the event, else discard event
  if((int)bjets.size() >= 4) 
  {
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


    // Check if it's a truth dijet
    if (abs(res_higgs[0].m() -120) > 0.1)
      res_weight = 0;
    if (abs(res_higgs[1].m() -120) > 0.1)
      res_weight = 0;

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

  } // /4jets


  // ************************************** Boosted Efficiencies ***********************************************

  double boost_weight = weightnorm;

  // Perform jet clustering with anti-kT
  // Note that here we use a small R clustering
  double const boostR=1.0; // To avoid overlapping b's as much as possible
  fastjet::JetDefinition boost_akt(fastjet::antikt_algorithm, boostR);

  // Get all the jets (no pt cut here)
  fastjet::ClusterSequence boost_cs_akt(fs, boost_akt);
  std::vector<fastjet::PseudoJet> fatjets = sorted_by_pt( boost_cs_akt.inclusive_jets()  );
  
  // We require 2 fatjets
  int const njet=2;
  if((int)fatjets.size() == njet) 
  {

    // Check if it's a truth higgs
    if (abs(fatjets[0].m() -120) > 0.1)
      boost_weight = 0;
    if (abs(fatjets[1].m() -120) > 0.1)
      boost_weight = 0;

    const fastjet::PseudoJet dihiggs= fatjets[0]+fatjets[1];

    FillHistogram("ptFatJet", boost_weight, fatjets[0].pt() );
    FillHistogram("ptFatJet", boost_weight, fatjets[1].pt() );

    FillHistogram("ptFatJet1", boost_weight, fatjets[0].pt() );
    FillHistogram("ptFatJet2", boost_weight, fatjets[1].pt() );

    FillHistogram("ptFatJetFatJet", boost_weight, fatjets[0].pt() + fatjets[1].pt());
  }
  


  // ************************************* MVA Output **********************************************************

  // Now save the ntuples to be used by the TMVA or the ANNs
  //   In the UCL analysis they use
  //   
  //   m, y, pT of the 4b system and masses of the two dijets
  //   3 decay angles (in resp. rest frames) & 2 angles between decay planes


  // This is for the UCL-like strategy
  // sabe mass, pt and y of th 4b system
  // the two dijet masses
  // and all independent angular distances between the four b jets
  // totalNTuple<<"# signal source m4b  pt4b y4b mHiggs1  mHiggs2 DeltaR_b1b2  DeltaR_b1b3  DeltaR_b1b4  DeltaR_b2b3  DeltaR_b2b4  DeltaR_b3b4 "<<std::endl;
  outputNTuple <<signal <<"\t"<<GetSample()<<"\t"<<event_weight<<std::endl; 
  // Other combinations of kinematical variables could also be useful
  // Need to investigate the kinematics of the 4b final state in more detail

  // Pass event
  Pass(event_weight);

}
