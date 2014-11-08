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
  // Higgs histograms
  BookHistogram(new YODA::Histo1D(20, 0, 500), "pthh");
  BookHistogram(new YODA::Histo1D(20, 0, 600), "pth");

  const std::string tupleSpec = "# signal pthh ";
  outputNTuple<<tupleSpec<<std::endl;
}

void DurhamAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
  double event_weight = weightnorm;

  // Fetch jets
  std::vector<fastjet::PseudoJet> higgs_candidates;
  JetCluster_Durham(fs, higgs_candidates, event_weight);

  // Fails basic cuts
  if(event_weight<1e-30) return;

  // Form dihiggs candidate
  fastjet::PseudoJet dihiggs= higgs_candidates[0]+higgs_candidates[1];

  // *******************************************************************************

  // Histograms before cuts

  // ptH histogram
  FillHistogram("pth", event_weight, higgs_candidates[0].pt() );
  FillHistogram("pth", event_weight, higgs_candidates[1].pt() );


  // Histograms for the pt of the HH system
  FillHistogram("pthh", event_weight, dihiggs.pt() );

  // ************* CUTS ************************************************************

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

  // ************************************* MVA Output **********************************************************


  outputNTuple << signal <<"\t"<<GetSample()<<"\t"<<dihiggs.pt()<<std::endl;

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

  // Now look for substructure in each of these two dijets using the BDRS
  // mass-drop tagger
  for (int i = 0; i < njet; i++) 
  {
    // first recluster again with the large-R but with Cambridge-Aachen
    fastjet::ClusterSequence cs_sub(jets_akt[i].constituents(), CA10);
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
    Cut("massdrop", event_weight);
    event_weight=0.0;
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

  Cut("btag", initial_weight - event_weight);
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

