// ucl.cc

#include "ucl.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"

#include "YODA/Histo1D.h"

UCLAnalysis::UCLAnalysis(std::string const& sampleName):
Analysis("ucl", sampleName)
{
	const double ptb_min=0;
	const double ptb_max=600;
	const int nbin_ptb=20;

	// Higgs histograms
	BookHistogram(new YODA::Histo1D(20, 0, 500), "pthh");
	BookHistogram(new YODA::Histo1D(20, 0, 600), "pth");

	// b Jet histograms
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb1");
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb2");
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb3");
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb4");

	const std::string tupleSpec = "# signal source m4b  pt4b y4b mHiggs1  mHiggs2 DeltaR_b1b2  DeltaR_b1b3  DeltaR_b1b4  DeltaR_b2b3  DeltaR_b2b4  DeltaR_b3b4";
	outputNTuple<<tupleSpec<<std::endl;
}

void UCLAnalysis::Analyse(bool const& signal, finalState const& fs)
{
	double event_weight = 0.0;

	// Fetch jets
	std::vector<fastjet::PseudoJet> bjets_unsort;
	JetCluster_UCL(fs, bjets_unsort, event_weight);

	// Fails cuts
	if(event_weight<1e-30) return;

	// First of all, after basic selection, require that all four b jets are above 40 GeV
	double const pt_bjet_ucl = 40.0;
	// they should also be in central rapodity, |eta| < 2.5
	double const eta_bjet_ucl = 2.5;

	// Fill the histograms for the pt of the b jets before 
	// the corresponding kinematical cuts
	vector<fastjet::PseudoJet> bjets = sorted_by_pt(bjets_unsort);
	FillHistogram("ptb1", event_weight, bjets.at(0).pt() );
	FillHistogram("ptb2", event_weight, bjets.at(1).pt() );
	FillHistogram("ptb3", event_weight, bjets.at(2).pt() );
	FillHistogram("ptb4", event_weight, bjets.at(3).pt() );

	int const njet=4;
	// Restrict to the four leading jets in the event
	for(unsigned ijet=0; ijet<njet;ijet++)
	{
		if(bjets.at(ijet).pt() < pt_bjet_ucl || 
			fabs( bjets.at(ijet).eta() ) > eta_bjet_ucl) 
			{
				Cut("bjet_kin", event_weight);	// Kinematics cut on b-jets
				return;
			}
	}

	// The next step is to apply the dijet selection
	// Get the pairing that minimizes |m_dj1 - m_dj2|
	double dijet_mass[njet][njet];
	for(unsigned ijet=0;ijet<njet;ijet++)
		for(unsigned jjet=0;jjet<njet;jjet++)
		{
			// Compute jet masses
			const fastjet::PseudoJet sum = bjets[ijet] + bjets[jjet];
			dijet_mass[ijet][jjet] = sum.m();
		}

	double mdj_diff_min = 1e20; // Some large number to begin
	unsigned jet1_id1=10,jet1_id2=10,jet2_id1=10,jet2_id2=10;

	for(unsigned ijet=0;ijet<njet;ijet++)
		for(unsigned jjet=ijet+1;jjet<njet;jjet++)
			for(unsigned ijet2=0;ijet2<njet;ijet2++)
				for(unsigned jjet2=ijet2+1;jjet2<njet;jjet2++)
				{
					double mdj1 = dijet_mass[ijet][jjet];
					double mdj2 = dijet_mass[ijet2][jjet2];
					double min_dj = fabs(mdj1 - mdj2);
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
	const fastjet::PseudoJet higgs1 = bjets.at( jet1_id1) + bjets.at( jet1_id2); 
	const fastjet::PseudoJet higgs2 = bjets.at( jet2_id1) + bjets.at( jet2_id2);

	// Histograms before cut in pt Higgs candidates
	FillHistogram("pth", event_weight, higgs1.pt() );
	FillHistogram("pth", event_weight, higgs2.pt() );


	// Now, the pt of the dijet Higgs candidates must be above 150 GeV
	const double pt_dijet_ucl=150.0;
	if( higgs1.pt() < pt_dijet_ucl || higgs2.pt() < pt_dijet_ucl ) // Was bugged, to higgs1 in both cases
	{
		Cut("pT_dijet", event_weight);	// Kinematics cut on b-jets 
		return;
	}

	// These two dijets cannot be too far in DeltaR
	// Check exactly the cut used: deltaR or delta_eta?
	const double delta_eta_dijet_ucl=1.5;
	const double delta_eta_dijet = fabs(higgs1.eta()- higgs2.eta());
	if(delta_eta_dijet > delta_eta_dijet_ucl) 
	{
		Cut("deltaR_dijet", event_weight);
		return;
	}

	// Higgs mass window condition
	const double mass_diff1 = fabs(higgs1.m()-m_higgs)/m_higgs;
	const double mass_diff2 = fabs(higgs2.m()-m_higgs)/m_higgs;
	if( mass_diff1 > mass_resolution || mass_diff2 > mass_resolution ) 
	{
		Cut("Higgs_massRes", event_weight);
		return;
	}
	// Histograms for the pt of the HH system
	// no cuts are applied on this variable
	const fastjet::PseudoJet dihiggs= higgs1+higgs2;
	FillHistogram("pthh", event_weight, dihiggs.pt() );

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
	outputNTuple <<signal <<"\t"<<GetSample()<<"\t"<<dihiggs.m()<<"\t"<<dihiggs.pt()<<"\t"<<dihiggs.rapidity()<<"\t"<<
	higgs1.m()<<"\t"<<higgs2.m()<<"\t"<<
	bjets.at(0).delta_R(bjets.at(1))<<"\t"<<
	bjets.at(0).delta_R(bjets.at(2))<<"\t"<<
	bjets.at(0).delta_R(bjets.at(3))<<"\t"<<
	bjets.at(1).delta_R(bjets.at(2))<<"\t"<<
	bjets.at(1).delta_R(bjets.at(2))<<"\t"<<
	bjets.at(2).delta_R(bjets.at(3))<<std::endl; 
	// Other combinations of kinematical variables could also be useful
	// Need to investigate the kinematics of the 4b final state in more detail

	// Pass event
	Pass(event_weight);

}

/*
This routine read the event kinematics and performs the jet clustering
It also checkes that energy momentum is conserved event by event
This applies for small R jet clustering with the anti-kt algorithm
 */

void UCLAnalysis::JetCluster_UCL(finalState const& particles, std::vector<fastjet::PseudoJet>& bjets, double& event_weight)
{
  // Perform jet clustering with anti-kT
  // Note that here we use a small R clustering
  static double const jetR=0.5; // To avoid overlapping b's as much as possible
  fastjet::JetDefinition akt(fastjet::antikt_algorithm, jetR);
  // Cluster all particles
  // The cluster sequence has to be saved to be used for jet substructure
  fastjet::ClusterSequence cs_akt(particles, akt);
  // Get all the jets (no pt cut here)
  std::vector<fastjet::PseudoJet> jets_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
  
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
  
  // Check energy-momentum conservation
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

  // Now Initialize the event weight
event_weight=1.0;

  // We require at least 4 jets in the event, else discard event
int const njet=4;
if(jets_akt.size() < njet) 
{
	Cut("4jet",event_weight);
	event_weight=0;
	return;
}

  // By looking at the jet constituents
  // we can simulate the effects of b tagging

  // Loop over the 4 hardest jets in event only
for(unsigned ijet=0; ijet<njet;ijet++)
	if( BTagging(jets_akt[ijet]) )   // Check if at least one of its constituents is a b quark
	{
		bjets.push_back(jets_akt.at(ijet));
		event_weight *= btag_prob; // Account for b tagging efficiency
	}
	else // Else, account for the fake b-tag probabililty
	{
		bjets.push_back(jets_akt.at(ijet));
		event_weight *= btag_mistag;
	}

	// cut from btagging
	Cut("btag", 1.0 - event_weight);

} 

// ----------------------------------------------------------------------------------

bool UCLAnalysis::BTagging( fastjet::PseudoJet const& jet ) const
{
	// Cuts for the b-jet candidates for b-tagging
	double const pt_btagging=15;

	// Get the jet constituents
	const std::vector<fastjet::PseudoJet>& jet_constituents = jet.constituents();

	// Loop over constituents and look for b quarks
	// also b quarks must be above some minimum pt
	for(size_t i=0; i<jet_constituents.size(); i++)
	{
		// Flavour of jet constituent
		const int userid= jet_constituents.at(i).user_index();
		const double pt_bcandidate = jet_constituents.at(i).pt();

		if(abs(userid) ==5 )
			if( pt_bcandidate > pt_btagging)
		  		return true;
	}

 	return false; // no b-jets found
}


