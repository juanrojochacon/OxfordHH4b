// oxford_boost_vr.cc

#include "oxford_boost_vr.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

#include "YODA/Histo1D.h"

using namespace fastjet::contrib;

// Set VR parameters
static double const jet_Rmax	=1.0;
static double const jet_Rmin	=0.1;
static double const jet_Rho	=500.;

//Instantiate VR plugin
static const VariableRPlugin lvjet_pluginAKT(jet_Rho, jet_Rmin, jet_Rmax, VariableRPlugin::AKTLIKE);
static const fastjet::JetDefinition VR_AKT(&lvjet_pluginAKT);


OxfordBoostVRAnalysis::OxfordBoostVRAnalysis(std::string const& sampleName):
Analysis("oxford_boost_vr", sampleName)
{
  // Plotting parameters
  const double ptfj_min=0;
  const double ptfj_max=900;

  const double DeltaRmin = 1;
  const double DeltaRmax = 5;
  
  // *********************** preCut **************************

  // Fat Jet histograms
  BookHistogram(new YODA::Histo1D(20, 0, 20), "nJets_preCut");
  
  BookHistogram(new YODA::Histo1D(20, ptfj_min, ptfj_max), "ptfj1_preCut");
  BookHistogram(new YODA::Histo1D(20, ptfj_min, ptfj_max), "ptfj2_preCut");

  BookHistogram(new YODA::Histo1D(20, 0, 200), "mfj1_preCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "mfj2_preCut");
  
  BookHistogram(new YODA::Histo1D(20, 0, 200), "split12_fj1_preCut");
  BookHistogram(new YODA::Histo1D(20, 0, 200), "split12_fj2_preCut");
  BookHistogram(new YODA::Histo1D(20, 0, 1), "tau21_fj1_preCut");
  BookHistogram(new YODA::Histo1D(20, 0, 1), "tau21_fj2_preCut");
  
  // 2 fat jet system histograms
  BookHistogram(new YODA::Histo1D(20, 200, 1500), "m2fj_preCut");
  BookHistogram(new YODA::Histo1D(20, -2.5, 2.5), "y2fj_preCut");
  BookHistogram(new YODA::Histo1D(20, 200, 1500), "pT2fj_preCut");

  BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_fj1fj2_preCut");

  // ************************* postCut ********************************

  // Fat Jet histograms
  BookHistogram(new YODA::Histo1D(20, 0, 20), "nJets_postCut");
  
  BookHistogram(new YODA::Histo1D(40, ptfj_min, ptfj_max), "ptfj1_postCut");
  BookHistogram(new YODA::Histo1D(40, ptfj_min, ptfj_max), "ptfj2_postCut");

  BookHistogram(new YODA::Histo1D(40, 0, 200), "mfj1_postCut");
  BookHistogram(new YODA::Histo1D(40, 0, 200), "mfj2_postCut");
  
  BookHistogram(new YODA::Histo1D(40, 0, 200), "split12_fj1_postCut");
  BookHistogram(new YODA::Histo1D(40, 0, 200), "split12_fj2_postCut");
  BookHistogram(new YODA::Histo1D(40, 0, 1), "tau21_fj1_postCut");
  BookHistogram(new YODA::Histo1D(40, 0, 1), "tau21_fj2_postCut");
  
  // 2 fat jet system histograms
  BookHistogram(new YODA::Histo1D(40, 200, 1500), "m2fj_postCut");
  BookHistogram(new YODA::Histo1D(40, -2.5, 2.5), "y2fj_postCut");
  BookHistogram(new YODA::Histo1D(40, 200, 1500), "pT2fj_postCut");

  BookHistogram(new YODA::Histo1D(40, DeltaRmin, DeltaRmax), "DeltaR_fj1fj2_postCut");

	// ************************* cutFlow/Ntuples ********************************

	const std::string tupleSpec = "# signal source m2fj pthh y2fj mHiggs1 mHiggs2 split12_Higgs1 split12_Higgs2 tau21_Higgs1 tau21_Higgs2 DeltaR_fj1fj2";
	outputNTuple<<tupleSpec<<std::endl;

	// Order cutflow
	Cut("Basic: Two fatjets", 0);
	Cut("Basic: Fatjet cuts", 1);
	Cut("Basic: 4 small-R jets", 2);
	Cut("Basic: 2 subjets for each fatjet", 3);
	Cut("Basic: bTagging", 4);
}

void OxfordBoostVRAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
	// Set initial weight
	double event_weight = weightnorm;

	// Perform jet clustering with VR anti-kT
	std::vector<fastjet::PseudoJet> fatjets;
	
	// Fetch jets and substructure information
	std::vector<double> split12_vec;
	std::vector<double> tau21_vec;
	JetCluster_LargeVR(fs, fatjets, split12_vec, tau21_vec, event_weight);


	// Fails cuts
	if(event_weight<1e-30) return;
	
	//------------------------------------------------------------------------------
	// pt ordering and basic cuts now done in the JetCluster_LargeVR(...) function
	//------------------------------------------------------------------------------

	// Construct the Higgs candidates
	const fastjet::PseudoJet higgs1 = fatjets.at(0);
	const fastjet::PseudoJet higgs2 = fatjets.at(1);

	// Histograms for the pt of the HH system
	// no cuts are applied on this variable
	const fastjet::PseudoJet dihiggs= higgs1+higgs2;

	int nJets = fatjets.size();
	// Record jet multiplicity before cuts
	FillHistogram("nJets_postCut", event_weight, nJets );

	FillHistogram("ptfj1_postCut", event_weight, fatjets.at(0).pt() );
	FillHistogram("ptfj2_postCut", event_weight, fatjets.at(1).pt() );

	FillHistogram("mfj1_postCut", event_weight, fatjets.at(0).m() );
	FillHistogram("mfj2_postCut", event_weight, fatjets.at(1).m() );

	FillHistogram("split12_fj1_postCut", event_weight, split12_vec.at(0) );
	FillHistogram("split12_fj2_postCut", event_weight, split12_vec.at(1) );
	
	FillHistogram("tau21_fj1_postCut", event_weight, tau21_vec.at(0) );
	FillHistogram("tau21_fj2_postCut", event_weight, tau21_vec.at(1) );
	
	FillHistogram("m2fj_postCut", event_weight, dihiggs.m() );
	FillHistogram("pT2fj_postCut", event_weight, dihiggs.pt() );
	FillHistogram("y2fj_postCut", event_weight, dihiggs.rapidity() );
	FillHistogram("DeltaR_fj1fj2_postCut", event_weight, fatjets.at(0).delta_R(fatjets.at(1)) );

	// Now save the ntuples to be used by the TMVA or the ANNs
	//"# signal source m2fj pthh y2fj mHiggs1 mHiggs2 DeltaR_fj1fj2"
	outputNTuple <<signal <<"\t"<<GetSample()<<"\t"<<dihiggs.m()<<"\t"<<dihiggs.pt()<<"\t"<<dihiggs.rapidity()<<"\t"<<
	higgs1.m()<<"\t"<<higgs2.m()<<"\t"<<
	split12_vec.at(0)<<"\t"<<split12_vec.at(1)<<"\t"<<
	tau21_vec.at(0)<<"\t"<<tau21_vec.at(1)<<"\t"<<
	fatjets.at(0).delta_R(fatjets.at(1))<<std::endl; 

	// Pass event
	Pass(event_weight);
}

/*
This routine read the event kinematics and performs the jet clustering
It also checkes that energy momentum is conserved event by event
This applies for small R jet clustering with the anti-kt algorithm
 */

void OxfordBoostVRAnalysis::JetCluster_LargeVR(finalState const& fs, std::vector<fastjet::PseudoJet>& fatjets, std::vector<double>& split12_vec, std::vector<double>& tau21_vec, double& event_weight)
{
 
  // Get all the jets (no pt cut here)
  fastjet::ClusterSequence cs_akt(fs, VR_AKT);
  std::vector<fastjet::PseudoJet> jets_vr_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
  VerifyFourMomentum(jets_vr_akt);
  
  int nJets = jets_vr_akt.size();
  // We require at least 2 fatjets in the event, else discard event
  int const njet=2;
  if(nJets < njet) 
  {
	  Cut("Basic: Two fatjets",event_weight);
	  event_weight=0;
	  return;
  }
  
  // Calculate some substructure variables
  split12_vec = SplittingScales( jets_vr_akt );
  tau21_vec = NSubjettiness( jets_vr_akt, jet_Rmax, jet_Rmin, jet_Rho );
  
  const fastjet::PseudoJet dihiggs= jets_vr_akt.at(0)+jets_vr_akt.at(1);
  
  // Record jet multiplicity before cuts
  FillHistogram("nJets_preCut", event_weight, nJets );

  // Fill the histograms for the pt of the fat jets before 
  // the corresponding kinematical cuts
  FillHistogram("ptfj1_preCut", event_weight, jets_vr_akt.at(0).pt() );
  FillHistogram("ptfj2_preCut", event_weight, jets_vr_akt.at(1).pt() );

  FillHistogram("mfj1_preCut", event_weight, jets_vr_akt.at(0).m() );
  FillHistogram("mfj2_preCut", event_weight, jets_vr_akt.at(1).m() );

  FillHistogram("split12_fj1_preCut", event_weight, split12_vec.at(0) );
  FillHistogram("split12_fj2_preCut", event_weight, split12_vec.at(1) );
  
  FillHistogram("tau21_fj1_preCut", event_weight, tau21_vec.at(0) );
  FillHistogram("tau21_fj2_preCut", event_weight, tau21_vec.at(1) );
  
  FillHistogram("m2fj_preCut", event_weight, dihiggs.m() );
  FillHistogram("pT2fj_preCut", event_weight, dihiggs.pt() );
  FillHistogram("y2fj_preCut", event_weight, dihiggs.rapidity() );
  FillHistogram("DeltaR_fj1fj2_preCut", event_weight, jets_vr_akt.at(0).delta_R(jets_vr_akt.at(1)) );
  
  
  //------------------------------------------------------------------------------------------------------
  // First of all, after basic selection, require that both jets are above 100 GeV
  double const pt_fatjet_ox = 100.0;
  // they should also be in central rapodity, |eta| < 2.5
  double const eta_fatjet_ox = 2.5;
  
  //Basic kinematic cuts on the jets
  int const nfatjet=2;
  // Restrict to the two leading jets in the event
  for(int ijet=0; ijet<nfatjet;ijet++)
  {
	  if(jets_vr_akt.at(ijet).pt() < pt_fatjet_ox || 
		  fabs( jets_vr_akt.at(ijet).eta() ) > eta_fatjet_ox) 
		  {
			  Cut("Basic: Fatjet cuts", event_weight);
			  event_weight=0;
			  return;
		  }
		  else fatjets.push_back(jets_vr_akt.at(ijet));
  }
  
  //------------------------------------------------------------------------------------------------------
  // Simulate subjets:
  // Perform small-R jet clustering with anti-kT
  static double const jetR=0.3; // To avoid overlapping b's as much as possible
  fastjet::JetDefinition jd_subjets(fastjet::antikt_algorithm, jetR);
  fastjet::ClusterSequence cs_subjets(fs, jd_subjets);

  std::vector<fastjet::PseudoJet> jets_akt = sorted_by_pt( cs_subjets.inclusive_jets()  );
  
  std::vector<fastjet::PseudoJet> subjets;
  
  //------------------------------------------------------------------------------------------------------
  //Consider only jets that pass basic kinematic cuts
  double const pt_subjet_ox = 25.0;
  // they should also be in central rapidity, |eta| < 2.5
  double const eta_subjet_ox = 2.5;
  
  //Basic kinematic cuts on the jets
  // Restrict to the four leading subjets in the event
  for(int ijet=0; ijet<(int)subjets.size();ijet++)
  {
	  if(subjets.at(ijet).pt() > pt_subjet_ox && 
		  fabs( jets_akt.at(ijet).eta() ) < eta_subjet_ox) 
		  {
			  subjets.push_back(jets_akt.at(ijet));
		  }
  }

  // We require at least 4 selected subjets in the event, else discard event
  int const nsubjet=4;
  if((int)subjets.size() < nsubjet) 
  {
	  Cut("Basic: 4 small-R jets",event_weight);
	  event_weight=0;
	  return;
  }

  std::vector<fastjet::PseudoJet> subjets_jet0;
  std::vector<fastjet::PseudoJet> subjets_jet1;
  
  get_assoc_trkjets( jets_vr_akt.at(0), jets_akt, subjets_jet0, false);
  get_assoc_trkjets( jets_vr_akt.at(1), jets_akt, subjets_jet1, false);
  
  // Require at least 2 subjets matched to each fatjet
  if((int)subjets_jet0.size() < 2 || (int)subjets_jet1.size() < 2) 
  {
	  Cut("Basic: 2 subjets for each fatjet",event_weight);
	  event_weight=0;
	  return;
  }
  
  std::vector<fastjet::PseudoJet> bjets_jet0;
  std::vector<fastjet::PseudoJet> bjets_jet1;
  
  // By looking at the jet constituents
  // we can simulate the effects of b tagging
  const double initial_weight = event_weight;
  for(unsigned int ijet=0; ijet<subjets_jet0.size(); ijet++){
  	if( BTagging(subjets_jet0[ijet]) )   // Check if at least two of its constituents are b quarks
  	{
  		bjets_jet0.push_back(subjets_jet0.at(ijet));
  		event_weight *= btag_prob; // Account for b tagging efficiency
  	}
  	else // Else, account for the fake b-tag probability
  	{
  		event_weight *= btag_mistag;
  	}
  }
  for(unsigned int ijet=0; ijet<subjets_jet1.size(); ijet++){
  	if( BTagging(subjets_jet1[ijet]) )   // Check if at least one of its constituents are b quarks
  	{
  		bjets_jet1.push_back(subjets_jet1.at(ijet));
  		double btag_prob = btag_eff( subjets_jet1.at(ijet).pt() );
  		event_weight *= btag_prob; // Account for b tagging efficiency
  	}
  	else if( CTagging(subjets_jet1[ijet])) {
		double ctag_prob = charm_eff( subjets_jet1.at(ijet).pt() );
		event_weight *= ctag_prob; // Account for c (mis-)tagging efficiency
  	}
  	else // Else, account for the fake b-tag probability
  	{	
		double mistag_prob = mistag_eff( subjets_jet1.at(ijet).pt() );
  		event_weight *= mistag_prob;
  	}
  }
  
  if( bjets_jet0.size() < 2 || bjets_jet1.size() < 2 ){
  
      Cut("Basic: bTagging", initial_weight - event_weight);
      event_weight=0;
      return;
  }
  
  return;
}


// ----------------------------------------------------------------------------------

bool OxfordBoostVRAnalysis::BTagging( fastjet::PseudoJet const& jet ) const
{
	// Cuts for the b-jet candidates for b-tagging
	double const pt_btagging=15.;

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
	
	if(countB>0) return true;

 	return false; // no b-jets found
}


bool OxfordBoostVRAnalysis::CTagging( fastjet::PseudoJet const& jet ) const
{
	// Cuts for the b-jet candidates for b-tagging
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

		if(abs(userid) ==5 )
			if( pt_ccandidate > pt_ctagging)
		  		countC++;
	}
	
	if(countC>0) return true;

 	return false; // no c-jets found
}
