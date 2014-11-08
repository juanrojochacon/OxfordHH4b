// durham.cc

#include "durham.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"

#include "YODA/Histo1D.h"

//-------------------------------------------------------------------
// Durham Analysis settings

static double const jetR_1p2=1.2; // To try to merge two b quarks into the same jet
static fastjet::JetDefinition CA10(fastjet::cambridge_algorithm, jetR_1p2);

// Set parameters of the mass drop tagger for jet substructure
// mu = 0.67 and y = 0.09 are the default choice in FastJet
static const double mu = 0.67;
static const double ycut = 0.09;
// The choice ycut = 0.15 is also recommended by Gavin to optimize S/sqrt(B)

//-------------------------------------------------------------------


DurhamAnalysis::DurhamAnalysis(std::string const& sampleName):
Analysis("durham", sampleName)
{
  // Plotting parameters
  const double ptfj_min=0;
  const double ptfj_max=900;
  const int nbin_ptfj=20;

  const double DeltaRmin = 0;
  const double DeltaRmax = 5;
  
  // *********************** preCut **************************

  // Fat Jet histograms
  BookHistogram(new YODA::Histo1D(nbin_ptfj, ptfj_min, ptfj_max), "ptfj1_preCut");
  BookHistogram(new YODA::Histo1D(nbin_ptfj, ptfj_min, ptfj_max), "ptfj2_preCut");

  BookHistogram(new YODA::Histo1D(20, 0, 200), "mfj1_preCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "mfj2_preCut");
  
  BookHistogram(new YODA::Histo1D(20, 0, 200), "split12_fj1_preCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "split12_fj2_preCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "tau21_fj1_preCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "tau21_fj2_preCut");
  
  // 2 fat jet system histograms
  BookHistogram(new YODA::Histo1D(20, 200, 1500), "m2fj_preCut");
  BookHistogram(new YODA::Histo1D(20, -2.5, 2.5), "y2fj_preCut");
  BookHistogram(new YODA::Histo1D(20, 200, 1500), "pT2fj_preCut");

  BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_fj1fj2_preCut");

  // ************************* postCut ********************************

  // Fat Jet histograms
  BookHistogram(new YODA::Histo1D(nbin_ptfj, ptfj_min, ptfj_max), "ptfj1_postCut");
  BookHistogram(new YODA::Histo1D(nbin_ptfj, ptfj_min, ptfj_max), "ptfj2_postCut");

  BookHistogram(new YODA::Histo1D(20, 0, 200), "mfj1_postCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "mfj2_postCut");
  
  BookHistogram(new YODA::Histo1D(20, 0, 200), "split12_fj1_postCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "split12_fj2_postCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "tau21_fj1_postCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "tau21_fj2_postCut");
  
  // 2 fat jet system histograms
  BookHistogram(new YODA::Histo1D(20, 200, 1500), "m2fj_postCut");
  BookHistogram(new YODA::Histo1D(20, -2.5, 2.5), "y2fj_postCut");
  BookHistogram(new YODA::Histo1D(20, 200, 1500), "pT2fj_postCut");

  BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_fj1fj2_postCut");

  // ************************* cutFlow/Ntuples ********************************

  // CutFlow
  Cut("Basic: Two FatJets", 0);
  Cut("Basic: bTagging", 0);

  const std::string tupleSpec = "# signal source m2fj pthh y2fj mHiggs1 mHiggs2 split12_Higgs1 split12_Higgs2 tau21_Higgs1 tau21_Higgs2 DeltaR_fj1fj2";
  outputNTuple<<tupleSpec<<std::endl;

}

void DurhamAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
  double event_weight = weightnorm;
  const int njet = 2; // limit to 2 hardest jets

  // Fetch jets
  std::vector<fastjet::PseudoJet> fatjets;
  JetCluster_Durham(fs, fatjets, event_weight);

  // Fails basic cuts
  if(event_weight<1e-30) return;

  // Form dihiggs candidate
  fastjet::PseudoJet fj2= fatjets[0]+fatjets[1];

   // Calculate some substructure variables
  std::vector<double> split12_vec = SplittingScales( fatjets );
  std::vector<double> tau21_vec = NSubjettiness( fatjets, jetR_1p2 );

  // *******************************************************************************

  // Histograms before cuts

  // Fill the histograms for the fat jets before the corresponding kinematical cuts
  FillHistogram("ptfj1_preCut", event_weight, fatjets.at(0).pt() );
  FillHistogram("ptfj2_preCut", event_weight, fatjets.at(1).pt() );

  FillHistogram("mfj1_preCut", event_weight, fatjets.at(0).m() );
  FillHistogram("mfj2_preCut", event_weight, fatjets.at(1).m() );

  FillHistogram("split12_fj1_preCut", event_weight, split12_vec.at(0) );
  FillHistogram("split12_fj2_preCut", event_weight, split12_vec.at(1) );
  
  FillHistogram("tau21_fj1_preCut", event_weight, tau21_vec.at(0) );
  FillHistogram("tau21_fj2_preCut", event_weight, tau21_vec.at(1) );

  FillHistogram("m2fj_preCut", event_weight, fj2.m() );
  FillHistogram("y2fj_preCut", event_weight, fj2.rapidity() );
  FillHistogram("pT2fj_preCut", event_weight, fj2.pt() );

  FillHistogram("DeltaR_fj1fj2_preCut", event_weight, fatjets[0].delta_R(fatjets[1]) );

  // ************* CUTS ************************************************************

  // Now look for substructure in each of these two dijets using the BDRS mass-drop tagger
  std::vector<fastjet::PseudoJet> higgs_candidates;
  for (int i = 0; i < njet; i++) 
  {
    // first recluster again with the large-R but with Cambridge-Aachen
    fastjet::ClusterSequence cs_sub(fatjets[i].constituents(), CA10);
    // get hardest jet
    fastjet::PseudoJet ca_jet = sorted_by_pt(cs_sub.inclusive_jets())[0];

    // now run mass drop tagger
    // parameters are specified in settings.h
    fastjet::MassDropTagger md_tagger(mu, ycut);
    fastjet::PseudoJet tagged_jet = md_tagger(ca_jet);

    // If tagging succesful - declare as Higgs candidate
    if (tagged_jet != 0 )  
      higgs_candidates.push_back(tagged_jet);
  }

  // If we don't have a mass-drop tag in each of the two leading large-R jets
  // discard the event
  if(higgs_candidates.size()!=2) 
  {
    Cut("mass-drop", event_weight);
    event_weight=0.0;
    return;
  }

  // Now the cut on the pt of the Higgs candidates
  double const pt_largeRjet = 200.0;
  for(size_t ijet=0; ijet<higgs_candidates.size(); ijet++)
    if(higgs_candidates[ijet].pt() < pt_largeRjet) 
    {
      Cut("jet_pT", event_weight);
      event_weight=0;
      return;
    }

  // Require that the Higgs candidates are in acceptance
  double const eta_bjet_durham = 2.5;
  for(size_t ijet=0; ijet<higgs_candidates.size(); ijet++)
    if(fabs( higgs_candidates.at(ijet).eta() ) > eta_bjet_durham) 
      return Cut("bjet_acceptance", event_weight);

  // Same as in the UCL analysis
  // require that these two leading fat jets are not too separated in rapidity
  const double delta_eta_dijet_fatjet=1.5;
  double delta_eta_dijet = fabs(higgs_candidates.at(0).eta()- higgs_candidates.at(1).eta());
  if(delta_eta_dijet > delta_eta_dijet_fatjet) return Cut("delta_eta_dijet", event_weight);

  // Higgs mass window condition
  const double mass_diff1 = fabs(higgs_candidates.at(0).m()-m_higgs)/m_higgs;
  const double mass_diff2 = fabs(higgs_candidates.at(1).m()-m_higgs)/m_higgs;
  if( mass_diff1 > mass_resolution || mass_diff2 > mass_resolution ) return  Cut("Higgs_window", event_weight);


  // ************************************ Post cut fills ***********************************************

  // Histograms after cuts

  // Fill the histograms for the fat jets before the corresponding kinematical cuts
  FillHistogram("ptfj1_postCut", event_weight, fatjets.at(0).pt() );
  FillHistogram("ptfj2_postCut", event_weight, fatjets.at(1).pt() );

  FillHistogram("mfj1_postCut", event_weight, fatjets.at(0).m() );
  FillHistogram("mfj2_postCut", event_weight, fatjets.at(1).m() );

  FillHistogram("split12_fj1_postCut", event_weight, split12_vec.at(0) );
  FillHistogram("split12_fj2_postCut", event_weight, split12_vec.at(1) );
  
  FillHistogram("tau21_fj1_postCut", event_weight, tau21_vec.at(0) );
  FillHistogram("tau21_fj2_postCut", event_weight, tau21_vec.at(1) );

  FillHistogram("m2fj_postCut", event_weight, fj2.m() );
  FillHistogram("y2fj_postCut", event_weight, fj2.rapidity() );
  FillHistogram("pT2fj_postCut", event_weight, fj2.pt() );

  FillHistogram("DeltaR_fj1fj2_postCut", event_weight, fatjets[0].delta_R(fatjets[1]) );

  // ************************************* MVA Output **********************************************************


  outputNTuple << signal <<"\t"<<GetSample()<<"\t"<<fj2.pt()<<std::endl;

  // Pass event
  Pass(event_weight);

}

/*
This routine read the event kinematics and performs the jet clustering
It also checkes that energy momentum is conserved event by event
This applies for small R jet clustering with the anti-kt algorithm
*/

void DurhamAnalysis::JetCluster_Durham(finalState const& particles, std::vector<fastjet::PseudoJet>& higgs_candidates, double& event_weight)
{  
  // Perform jet clustering with anti-kT
  // jetR is defined in settings.h
  // Here we use a large-R clustering, R=1.2
  fastjet::JetDefinition akt(fastjet::antikt_algorithm, jetR_1p2);
  fastjet::ClusterSequence cs_akt(particles, akt);

  std::vector<fastjet::PseudoJet> jets_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
  VerifyFourMomentum(jets_akt);

  // We require at least 2 large-R jets in the event, else discard event
  int const njet=2;
  if((int)jets_akt.size() < njet) 
  {
    Cut("Basic: Two FatJets", event_weight);
    event_weight=0;
    return;
  }

  // Get the jet constituents
  int const nb_fatjet=2;
  const double initial_weight = event_weight;
  for (int i = 0; i < njet; i++) 
  {
    // Tag b's
    const int nb = BTagging(jets_akt[i]);

    // Now assign the event weight
    // This can be refined, for example requiring for fakes to have subjets with
    // pt above some threshold
    if(nb >= nb_fatjet) event_weight *= pow(btag_prob,2.0);
      if(nb ==1 ) event_weight *= (btag_prob*btag_mistag);
        if(nb ==0 ) event_weight *= pow(btag_mistag,2.0);
  }

  Cut("Basic: bTagging", initial_weight - event_weight);
  return;

} 

// ----------------------------------------------------------------------------------

// Now we need to implement b-tagging
// Here different options possible
// To begin with, most agressive option
// Require only two b quarks in the fat jet
// each with pt > 40 GeV
// No restruction on how close they can be in angle
int DurhamAnalysis::BTagging(fastjet::PseudoJet const& jet) const
{
  // pT cut in Durham approach
  double const pt_btagging_largeR=40.0;

  const std::vector<fastjet::PseudoJet>& jet_constituents = jet.constituents();

  // Loop over constituents and look for b quarks
  // These b quarks must be hard enough
  int nb=0;
  for(size_t i=0;i<jet_constituents.size(); i++)
  {
    const int userid= jet_constituents.at(i).user_index();
    const double pt_bcandidate = jet_constituents.at(i).pt();

    if(abs(userid) ==5 )
      if( pt_bcandidate > pt_btagging_largeR) // Cut in pt
        nb++;
  }

  return nb;
}

