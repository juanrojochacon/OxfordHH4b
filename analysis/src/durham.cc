// ucl.cc

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

//-------------------------------------------------------------------

	
DurhamAnalysis::DurhamAnalysis():
Analysis("durham")
{
	// Higgs histograms
	BookHistogram(new YODA::Histo1D(20, 0, 500), "pthh");
	BookHistogram(new YODA::Histo1D(20, 0, 600), "pth");
	
 	tupleSpec = "# signal source pthh ";
 	totalNTuple<<tupleSpec<<std::endl;
}

void DurhamAnalysis::Analyse(string const& sampleID, bool const& signal, finalState const& fs)
{
	double event_weight = 0.0;

	// Fetch jets
	std::vector<fastjet::PseudoJet> higgs_candidates;
	JetCluster_Durham(fs, sampleID, higgs_candidates, event_weight);

	// Fails cuts
	if(event_weight<1e-30) return;

	// Require that the Higgs candidates are in acceptance
	double const eta_bjet_durham = 2.5;

	int const njet=2;
	for(unsigned ijet=0; ijet<njet;ijet++){
	if(fabs( higgs_candidates.at(ijet).eta() ) > eta_bjet_durham) return;
	}

	// Same as in the UCL analysis
	// require that these two leading fat jets are not too separated in rapidity
	const double delta_eta_dijet_fatjet=1.5;
	double delta_eta_dijet = fabs(higgs_candidates.at(0).eta()- higgs_candidates.at(1).eta());
	if(delta_eta_dijet > delta_eta_dijet_fatjet) return;

	// Higgs mass window condition
	const double mass_diff1 = fabs(higgs_candidates.at(0).m()-m_higgs)/m_higgs;
	const double mass_diff2 = fabs(higgs_candidates.at(1).m()-m_higgs)/m_higgs;
	if( mass_diff1 > mass_resolution || mass_diff2 > mass_resolution ) return;

	// Histograms for the pt of the HH system
	// no cuts are applied on this variable
	fastjet::PseudoJet dihiggs= higgs_candidates.at(0)+higgs_candidates.at(1);
	FillHistogram("pthh", sampleID, event_weight, dihiggs.pt() );

	totalNTuple << signal <<"\t" << sampleID <<"\t"<<dihiggs.pt()<<endl;
  	sampleNTuple << signal <<"\t" << sampleID <<"\t"<<dihiggs.pt()<<endl;

    // Increment passed counter
	nPassed++;
	passedWeight += event_weight;
	
}

/*
This routine read the event kinematics and performs the jet clustering
It also checkes that energy momentum is conserved event by event
This applies for small R jet clustering with the anti-kt algorithm
 */

void DurhamAnalysis::JetCluster_Durham(finalState const& particles, string const& sampleID, std::vector<fastjet::PseudoJet>& higgs_candidates, double& event_weight)
{  
  // Perform jet clustering with anti-kT
  // jetR is defined in settings.h
  // Here we use a large-R clustering, R=1.2
  fastjet::JetDefinition akt(fastjet::antikt_algorithm, jetR_1p2);
  // Cluster all particles
  // The cluster sequence has to be saved to be used for jet substructure
  fastjet::ClusterSequence cs_akt(particles, akt);
  // Get all the jets (no pt cut here)
  vector<fastjet::PseudoJet> jets_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
  
  // Check again four-momentum conservation, this time applied to jets
  // formed from the clustering of quarks and gluons (and beam remnants as well)
  double px_tot=0;
  double py_tot=0;
  double pz_tot=0;
  double E_tot=0;
  for(unsigned ij=0;ij<jets_akt.size();ij++){
    px_tot+= jets_akt.at(ij).px();
    py_tot+= jets_akt.at(ij).py();
    pz_tot+= jets_akt.at(ij).pz();
    E_tot+= jets_akt.at(ij).E();
  }
  if( fabs(px_tot) > tol_emom || fabs(py_tot)  > tol_emom 
      || fabs(pz_tot)  > tol_emom || fabs(E_tot-Eref)  > tol_emom ){
    std::cout<<"\n ********************************************************************** \n"<<std::endl;
    std::cout<<"No conservation of energy in Pythia after shower and jet reconstruction "<<std::endl;
    std::cout<<"px_tot = "<<px_tot<<std::endl;
    std::cout<<"py_tot = "<<py_tot<<std::endl;
    std::cout<<"pz_tot = "<<pz_tot<<std::endl;
    std::cout<<"E_tot, Eref = "<<E_tot<<" "<<Eref<<std::endl;
    exit(-10);
    std::cout<<"\n ********************************************************************** \n"<<std::endl;
  }
  
  // Now initialize the event weight
  event_weight=1.0;

  // We require at least 2 large-R jets in the event, else discard event
  int const njet=2;
  if(jets_akt.size() < njet) {
    event_weight=0;
    return;
  }

  // The two hardest jets must have pt > 200 GeV
  // Before doing the cut fill the ptH histogram
  FillHistogram("pth", sampleID, event_weight, jets_akt.at(0).pt() );
  FillHistogram("pth", sampleID, event_weight, jets_akt.at(1).pt() );

  // Now the cut on the pt of the Higgs candidates
  double const pt_largeRjet = 200.0;
  for(unsigned ijet=0;ijet<njet;ijet++){
    if(jets_akt.at(ijet).pt() < pt_largeRjet) {
      event_weight=0;
      return;
    }
  }

  // Now look for substructure in each of these two dijets using the BDRS
  // mass-drop tagger
  
  for (unsigned i = 0; i < njet; i++) {
    
    // first recluster again with the large-R but with Cambridge-Aachen
    fastjet::ClusterSequence cs_sub(jets_akt[i].constituents(), CA10);
    // get hardest jet
    fastjet::PseudoJet ca_jet = sorted_by_pt(cs_sub.inclusive_jets())[0];
    
    // now run mass drop tagger
    // parameters are specified in settings.h
    fastjet::MassDropTagger md_tagger(mu, ycut);
    fastjet::PseudoJet tagged_jet = md_tagger(ca_jet);
    
    // If tagging succesful - declare as Higgs candidate
    if (tagged_jet != 0 )  higgs_candidates.push_back(tagged_jet);
    
  }

  // If we don't have a mass-drop tag in each of the two leading large-R jets
  // discard the event
  if(higgs_candidates.size()!=2) {
    event_weight=0.0;
    return;
  }

  // Now we need to implement b-tagging
  // Here different options possible
  // To begin with, most agressive option
  // Require only two b quarks in the fat jet
  // each with pt > 40 GeV
  // No restruction on how close they can be in angle
  double const pt_btagging_largeR=40.0;

  // Get the jet constituents
  int const nb_fatjet=2;
  for (unsigned i = 0; i < njet; i++) {

    vector<fastjet::PseudoJet> jet_constituents = jets_akt.at(i).constituents();

    // number of b quarks
    int nb=0;

    // Loop over constituents and look for b quarks
    // These b quarks must be hard enough
    for(unsigned i=0;i<jet_constituents.size(); i++){
      int userid= jet_constituents.at(i).user_index();
      double pt_bcandidate = jet_constituents.at(i).pt();
      if(abs(userid) ==5 ){
	// Cut in pt
	if( pt_bcandidate > pt_btagging_largeR){
	  nb++;
	}
      }
    }
    // Now assign the event weight
    // This can be refined, for example requiring for fakes to have subjets with
    // pt above some threshold
    if(nb >= nb_fatjet) event_weight *= pow(btag_prob,2.0);
    if(nb ==1 ) event_weight *= (btag_prob*btag_mistag);
    if(nb ==0 ) event_weight *= pow(btag_mistag,2.0);
  }

  //std::cout<<"jet_clustering_analysis_largeR completed"<<std::endl;

  return;

} 

// ----------------------------------------------------------------------------------

