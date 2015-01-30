// oxford_res_fr.cc

#include "oxford_res_fr.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

#include "YODA/Histo1D.h"

using namespace fastjet::contrib;

OxfordResFRAnalysis::OxfordResFRAnalysis(std::string const& sampleName):
Analysis("oxford_res_fr", sampleName)
{
// ********************* Histogram settings******************

	const double DeltaRmin = 0;
	const double DeltaRmax = 5;

	const double DeltaPhimin = -3;
	const double DeltaPhimax = 3;

	const double DeltaEtamin = -2.5;
	const double DeltaEtamax = 2.5;

	const double ptb_min=0;
	const double ptb_max=600;
	const int nbin_ptb=20;

	// ********************* Histograms before cuts ***************************

	// Higgs histograms
	BookHistogram(new YODA::Histo1D(20, 0, 600), "ptdijet_preCut");
	BookHistogram(new YODA::Histo1D(20, 0, 600), "ptdijet1_preCut");
	BookHistogram(new YODA::Histo1D(20, 0, 600), "ptdijet2_preCut");
	BookHistogram(new YODA::Histo1D(20, 0, 500), "pt4b_preCut");

	// Histograms of dijet systems
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_dijet_preCut");
	BookHistogram(new YODA::Histo1D(20, DeltaPhimin, DeltaPhimax), "DeltaPhi_dijet_preCut");
	BookHistogram(new YODA::Histo1D(20, DeltaEtamin, DeltaEtamax), "DeltaEta_dijet_preCut");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "mDijet_preCut");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "mDijet1_preCut");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "mDijet2_preCut");

	// b Jet histograms
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb1_preCut");
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb2_preCut");
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb3_preCut");
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb4_preCut");

	BookHistogram(new YODA::Histo1D(5, 0, 5), "truth_NbJets");
	BookHistogram(new YODA::Histo1D(5, 0, 5), "truth_NbConstituents");
	BookHistogram(new YODA::Histo1D(5, 0, 5), "truth_NcJets");
	BookHistogram(new YODA::Histo1D(5, 0, 5), "truth_NcConstituents");
	BookHistogram(new YODA::Histo1D(20, 0, 100), "truth_cPT");
	BookHistogram(new YODA::Histo1D(20, 0, 100), "truth_bPT");


	// 4b system histograms
	BookHistogram(new YODA::Histo1D(20, 200, 1200), "m4b_preCut");
	BookHistogram(new YODA::Histo1D(20, -2.5, 2.5), "y4b_preCut");

	// Dijet distance
	BookHistogram(new YODA::Histo1D(20, -3, 3), "dijet_dijet_deltaPhi_preCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "minBJetDeltaR_preCut");

	// Delta-b histograms
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b1b2_preCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b1b3_preCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b1b4_preCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b2b3_preCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b2b4_preCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b3b4_preCut");


	// ********************* Histograms after ***************************

	// Higgs histograms
	BookHistogram(new YODA::Histo1D(20, 0, 600), "ptdijet_postCut");
	BookHistogram(new YODA::Histo1D(20, 0, 600), "ptdijet1_postCut");
	BookHistogram(new YODA::Histo1D(20, 0, 600), "ptdijet2_postCut");
	BookHistogram(new YODA::Histo1D(20, 0, 500), "pt4b_postCut");

	// Histograms of dijet systems
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_dijet_postCut");
	BookHistogram(new YODA::Histo1D(20, DeltaPhimin, DeltaPhimax), "DeltaPhi_dijet_postCut");
	BookHistogram(new YODA::Histo1D(20, DeltaEtamin, DeltaEtamax), "DeltaEta_dijet_postCut");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "mDijet_postCut");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "mDijet1_postCut");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "mDijet2_postCut");

	// b Jet histograms
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb1_postCut");
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb2_postCut");
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb3_postCut");
	BookHistogram(new YODA::Histo1D(nbin_ptb, ptb_min, ptb_max), "ptb4_postCut");

	// 4b system histograms
	BookHistogram(new YODA::Histo1D(20, 200, 1200), "m4b_postCut");
	BookHistogram(new YODA::Histo1D(20, -2.5, 2.5), "y4b_postCut");

	// Dijet distance
	BookHistogram(new YODA::Histo1D(20, -3, 3), "dijet_dijet_deltaPhi_postCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "minBJetDeltaR_postCut");

	// Delta-b histograms
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b1b2_postCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b1b3_postCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b1b4_postCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b2b3_postCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b2b4_postCut");
	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_b3b4_postCut");

	// ********************************************************************


	const std::string tupleSpec = "# signal source weight m4b  pt4b y4b mHiggs1  mHiggs2 DeltaR_b1b2  DeltaR_b1b3  DeltaR_b1b4  DeltaR_b2b3  DeltaR_b2b4  DeltaR_b3b4";
	outputNTuple<<tupleSpec<<std::endl;

	// Order cutflow
	Cut("Basic: Two dijets", 0);
	Cut("Basic: bTagging", 0);
	Cut("bJet pT/Eta", 0);
	Cut("diJet pT", 0);
	Cut("diJet DeltaR", 0);
	Cut("Higgs window", 0);
}

void OxfordResFRAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
	Analysis::Analyse(signal, weightnorm, fs);

	// Set initial weight
	double event_weight = weightnorm;

	// Fetch jets
	std::vector<fastjet::PseudoJet> bjets_unsort;
	JetCluster_SmallFR(fs, bjets_unsort, event_weight);

	// Fails cuts
	if(event_weight<1e-30) return;

	// Fill the histograms for the pt of the b jets before 
	// the corresponding kinematical cuts
	std::vector<fastjet::PseudoJet> bjets = sorted_by_pt(bjets_unsort);
	int const njet=4; 	// Restrict analysis to the four leading jets in the event

	// The next step is to apply the dijet selection
	// Get the pairing that minimizes |m_dj1 - m_dj2|
	double dijet_mass[njet][njet];
	for(int ijet=0;ijet<njet;ijet++)
		for(int jjet=0;jjet<njet;jjet++)
		{
			// Compute jet masses
			const fastjet::PseudoJet sum = bjets[ijet] + bjets[jjet];
			dijet_mass[ijet][jjet] = sum.m();
		}

	double mdj_diff_min = 1e20; // Some large number to begin
	int jet1_id1=10,jet1_id2=10,jet2_id1=10,jet2_id2=10;

	for(int ijet=0;ijet<njet;ijet++)
		for(int jjet=ijet+1;jjet<njet;jjet++)
			for(int ijet2=0;ijet2<njet;ijet2++)
				for(int jjet2=ijet2+1;jjet2<njet;jjet2++)
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
	std::vector<fastjet::PseudoJet> higgs;
	higgs.push_back(bjets.at( jet1_id1) + bjets.at( jet1_id2)); 
	higgs.push_back(bjets.at( jet2_id1) + bjets.at( jet2_id2));
	const fastjet::PseudoJet dihiggs= higgs[0]+higgs[1];
	// Sort by higgs pt
	higgs = sorted_by_pt(higgs);

	// Fetch minimum deltaR between b jets
	const double minDR11 = bjets[jet1_id1].delta_R(bjets[jet2_id1]);
	const double minDR12 = bjets[jet1_id1].delta_R(bjets[jet2_id2]);
	const double minDR21 = bjets[jet1_id2].delta_R(bjets[jet2_id1]);
	const double minDR22 = bjets[jet1_id2].delta_R(bjets[jet2_id2]);
	const double minDR = min(minDR11,min(minDR12,min(minDR21,minDR22)));

	// *******************************************************************************

	// Histograms before cuts

	FillHistogram("ptdijet_preCut", event_weight, higgs[0].pt() );
	FillHistogram("ptdijet_preCut", event_weight, higgs[1].pt() );

	FillHistogram("ptdijet1_preCut", event_weight, higgs[0].pt() );
	FillHistogram("ptdijet2_preCut", event_weight, higgs[1].pt() );

	FillHistogram("pt4b_preCut", event_weight, dihiggs.pt() );

	// dijet higgs candidates: delta<obs> between jets
	FillHistogram("DeltaR_dijet_preCut", event_weight, bjets.at(jet1_id1).delta_R(bjets.at(jet1_id2)) );
	FillHistogram("DeltaR_dijet_preCut", event_weight, bjets.at(jet2_id1).delta_R(bjets.at(jet2_id2)) );
	FillHistogram("DeltaPhi_dijet_preCut", event_weight, bjets.at(jet1_id1).delta_phi_to(bjets.at(jet1_id2)) );
	FillHistogram("DeltaPhi_dijet_preCut", event_weight, bjets.at(jet2_id1).delta_phi_to(bjets.at(jet2_id2)) );
	FillHistogram("DeltaEta_dijet_preCut", event_weight, bjets.at(jet1_id1).eta() - bjets.at(jet1_id2).eta() );
	FillHistogram("DeltaEta_dijet_preCut", event_weight, bjets.at(jet2_id1).eta() - bjets.at(jet2_id2).eta() );

	// Higgs candidate mass
	FillHistogram("mDijet_preCut", event_weight, higgs[0].m() );
	FillHistogram("mDijet_preCut", event_weight, higgs[1].m() );
	FillHistogram("mDijet1_preCut", event_weight, higgs[0].m() );
	FillHistogram("mDijet2_preCut", event_weight, higgs[1].m() );

	// bJet pT
	FillHistogram("ptb1_preCut", event_weight, bjets.at(0).pt() );
	FillHistogram("ptb2_preCut", event_weight, bjets.at(1).pt() );
	FillHistogram("ptb3_preCut", event_weight, bjets.at(2).pt() );
	FillHistogram("ptb4_preCut", event_weight, bjets.at(3).pt() );

	// dihiggs system
	FillHistogram("m4b_preCut", event_weight, dihiggs.m() );
	FillHistogram("y4b_preCut", event_weight, dihiggs.rapidity() );

	// Dijet distance
	FillHistogram("dijet_dijet_deltaPhi_preCut", event_weight, higgs[0].delta_phi_to(higgs[1]) );

	// Closest distance
	FillHistogram("minBJetDeltaR_preCut", event_weight, minDR );

	// bJet delta-R
	FillHistogram("DeltaR_b1b2_preCut", event_weight, bjets.at(0).delta_R(bjets.at(1)) );
	FillHistogram("DeltaR_b1b3_preCut", event_weight, bjets.at(0).delta_R(bjets.at(2)) );
	FillHistogram("DeltaR_b1b4_preCut", event_weight, bjets.at(0).delta_R(bjets.at(3)) );
	FillHistogram("DeltaR_b2b3_preCut", event_weight, bjets.at(1).delta_R(bjets.at(2)) );
	FillHistogram("DeltaR_b2b4_preCut", event_weight, bjets.at(1).delta_R(bjets.at(3)) );
	FillHistogram("DeltaR_b3b4_preCut", event_weight, bjets.at(2).delta_R(bjets.at(3)) );
	
	
// ************* CUTS ************************************************************

	// First of all, after basic selection, require that all four jets are above 25 GeV
	double const pt_bjet_ox = 25.0;
	// they should also be in central rapidity, |eta| < 2.5
	double const eta_bjet_ox = 2.5;

	// Perform cuts
	for(int ijet=0; ijet<njet;ijet++)
	{
		if(bjets.at(ijet).pt() < pt_bjet_ox || 
			fabs( bjets.at(ijet).eta() ) > eta_bjet_ox) 
			{
				Cut("Jet pT/Eta", event_weight);	// Kinematics cut on b-jets
				event_weight=0;
				return;
			}
	}

// *************************** Post cut fills **************************************

	FillHistogram("ptdijet_postCut", event_weight, higgs[0].pt() );
	FillHistogram("ptdijet_postCut", event_weight, higgs[1].pt() );

	FillHistogram("ptdijet1_postCut", event_weight, higgs[0].pt() );
	FillHistogram("ptdijet2_postCut", event_weight, higgs[1].pt() );

	FillHistogram("pt4b_postCut", event_weight, dihiggs.pt() );

	// dijet higgs candidates: delta<obs> between jets
	FillHistogram("DeltaR_dijet_postCut", event_weight, bjets.at(jet1_id1).delta_R(bjets.at(jet1_id2)) );
	FillHistogram("DeltaR_dijet_postCut", event_weight, bjets.at(jet2_id1).delta_R(bjets.at(jet2_id2)) );
	FillHistogram("DeltaPhi_dijet_postCut", event_weight, bjets.at(jet1_id1).delta_phi_to(bjets.at(jet1_id2)) );
	FillHistogram("DeltaPhi_dijet_postCut", event_weight, bjets.at(jet2_id1).delta_phi_to(bjets.at(jet2_id2)) );
	FillHistogram("DeltaEta_dijet_postCut", event_weight, bjets.at(jet1_id1).eta() - bjets.at(jet1_id2).eta() );
	FillHistogram("DeltaEta_dijet_postCut", event_weight, bjets.at(jet2_id1).eta() - bjets.at(jet2_id2).eta() );

	// Higgs candidate mass
	FillHistogram("mDijet_postCut", event_weight, higgs[0].m() );
	FillHistogram("mDijet_postCut", event_weight, higgs[1].m() );
	FillHistogram("mDijet1_postCut", event_weight, higgs[0].m() );
	FillHistogram("mDijet2_postCut", event_weight, higgs[1].m() );

	// bJet pT
	FillHistogram("ptb1_postCut", event_weight, bjets.at(0).pt() );
	FillHistogram("ptb2_postCut", event_weight, bjets.at(1).pt() );
	FillHistogram("ptb3_postCut", event_weight, bjets.at(2).pt() );
	FillHistogram("ptb4_postCut", event_weight, bjets.at(3).pt() );

	// dihiggs system
	FillHistogram("m4b_postCut", event_weight, dihiggs.m() );
	FillHistogram("y4b_postCut", event_weight, dihiggs.rapidity() );

	// Dijet distance
	FillHistogram("dijet_dijet_deltaPhi_postCut", event_weight, higgs[0].delta_phi_to(higgs[1]) );

	// Closest distance
	FillHistogram("minBJetDeltaR_postCut", event_weight, minDR );

	// bJet delta-R
	FillHistogram("DeltaR_b1b2_postCut", event_weight, bjets.at(0).delta_R(bjets.at(1)) );
	FillHistogram("DeltaR_b1b3_postCut", event_weight, bjets.at(0).delta_R(bjets.at(2)) );
	FillHistogram("DeltaR_b1b4_postCut", event_weight, bjets.at(0).delta_R(bjets.at(3)) );
	FillHistogram("DeltaR_b2b3_postCut", event_weight, bjets.at(1).delta_R(bjets.at(2)) );
	FillHistogram("DeltaR_b2b4_postCut", event_weight, bjets.at(1).delta_R(bjets.at(3)) );
	FillHistogram("DeltaR_b3b4_postCut", event_weight, bjets.at(2).delta_R(bjets.at(3)) );

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
	outputNTuple <<signal <<"\t"<<GetSample()<<"\t"<<event_weight<<"\t"<<dihiggs.m()<<"\t"<<dihiggs.pt()<<"\t"<<dihiggs.rapidity()<<"\t"<<
	higgs[0].m()<<"\t"<<higgs[1].m()<<"\t"<<
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

void OxfordResFRAnalysis::JetCluster_SmallFR(finalState const& particles, std::vector<fastjet::PseudoJet>& bjets, double& event_weight)
{
	
	static double const jetR=0.4; // To avoid overlapping b's as much as possible
	fastjet::JetDefinition akt(fastjet::antikt_algorithm, jetR);

	// Cluster all particles
	// The cluster sequence has to be saved to be used for jet substructure
	fastjet::ClusterSequence cs_akt(particles, akt);
	// Get all the jets (no pt cut here)
	std::vector<fastjet::PseudoJet> jets_fr_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
	VerifyFourMomentum(jets_fr_akt);

	// We require at least 4 jets in the event, else discard event
	int const njet=4;
	if((int)jets_fr_akt.size() < njet) 
	{
		Cut("Basic: Two dijets",event_weight);
		event_weight=0;
		return;
	}
	
	// By looking at the jet constituents
	// we can simulate the effects of b tagging

	int Nbjets=0;
	int Ncjets=0;

	// Loop over the 4 hardest jets in event only
	const double initial_weight = event_weight;
	for(int ijet=0; ijet<njet;ijet++) {
	  int bQuarks = BTagging(jets_fr_akt[ijet]);
	  int cQuarks = CTagging(jets_fr_akt[ijet]);
	  //int cQuarks = 0;
	  FillHistogram("truth_NbConstituents", 1, bQuarks+0.5 );
	  FillHistogram("truth_NcConstituents", 1, cQuarks+0.5 );
	  
	  if( bQuarks > 0 )   // Check if at least one of its constituents is a b quark
	    {
	      bjets.push_back(jets_fr_akt.at(ijet));
	      event_weight *= btag_prob; // Account for b tagging efficiency
	      Nbjets++;
	    }
	  else if ( cQuarks > 0 ) {
	    bjets.push_back(jets_fr_akt.at(ijet));
	    event_weight *= ctag_prob; // Include c-mis tag rate
	    Ncjets++;
	  }
	
	  else // Else, account for the fake b-tag probabililty
	    {
	      bjets.push_back(jets_fr_akt.at(ijet));
	      event_weight *= btag_mistag;
	    }
	  
	  FillHistogram("truth_NbJets", 1, Nbjets+0.5 );
	  FillHistogram("truth_NcJets", 1, Ncjets+0.5 );

	}

	// cut from btagging
	Cut("Basic: bTagging", initial_weight - event_weight);
} 

// ----------------------------------------------------------------------------------

int OxfordResFRAnalysis::BTagging( fastjet::PseudoJet const& jet ) 
{
	// Cuts for the b-jet candidates for b-tagging
	double const pt_btagging=15;

	// Get the jet constituents
	const std::vector<fastjet::PseudoJet>& jet_constituents = jet.constituents();
	int countB = 0;

	// Loop over constituents and look for b quarks
	// also b quarks must be above some minimum pt
	for(size_t i=0; i<jet_constituents.size(); i++)
	{
		// Flavour of jet constituent
		const int userid= jet_constituents.at(i).user_index();
		const double pt_bcandidate = jet_constituents.at(i).pt();

		if(abs(userid) ==5 )
		  {     FillHistogram("truth_bPT", 1, pt_bcandidate);
			if( pt_bcandidate > pt_btagging)
		  		countB++;
		  }
	}

 	return countB; // no b-jets found
}

//----------------------------------------------------------------------------------

int OxfordResFRAnalysis::CTagging( fastjet::PseudoJet const& jet ) 
{
  // Cuts for the c-quark candidates for c-(mis)tagging
  double const pt_ctagging=15.;
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
	{ FillHistogram("truth_cPT", 1, pt_ccandidate);
	  if( pt_ccandidate > pt_ctagging)
	    countC++;
	}
    }
  return countC;
}
