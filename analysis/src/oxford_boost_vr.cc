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
static double const jet_Rho	=250.;

// Set parameters for mass-drop tagger
// mu = 0.67 and y = 0.09 are the default choice in FastJet
static double const mu = 0.67;
static double const ycut = 0.09;
// The choice ycut = 0.15 is also recommended by Gavin to optimize S/sqrt(B)

//Instantiate VR plugin
static const VariableRPlugin lvjet_pluginAKT(jet_Rho, jet_Rmin, jet_Rmax, VariableRPlugin::AKTLIKE);
static const fastjet::JetDefinition VR_AKT(&lvjet_pluginAKT);

// CA for reclustering
double const R_CA_boost = 1.2;
static const fastjet::JetDefinition CA12(fastjet::cambridge_algorithm, R_CA_boost);


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
	Cut("Two dijets", 0);
	Cut("Higgs window", 0);
}

void OxfordBoostVRAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
	// Set initial weight
	double event_weight = weightnorm;

	// Fetch jets and substructure information
	std::vector<fastjet::PseudoJet> fatjets;
	fastjet::ClusterSequence cs_akt(fs, VR_AKT);

	JetCluster_LargeVR(cs_akt, fatjets, event_weight);

  //-----------------------------------------------
  // Mass drop tagger
  //-----------------------------------------------
  vector<fastjet::PseudoJet> higgs_candidates;
  for (unsigned i = 0; i < fatjets.size(); i++) {
  
    // first recluster with some large CA (needed for mass-drop)
    // Would it help to switch to CA VR to avoid this step?
    fastjet::ClusterSequence cs_sub(fatjets[i].constituents(), CA12);
    
    // next get hardest jet
    fastjet::PseudoJet ca_jet = sorted_by_pt(cs_sub.inclusive_jets())[0];
    
    // now run mass drop tagger
    fastjet::MassDropTagger md_tagger(mu, ycut);
    fastjet::PseudoJet tagged_jet = md_tagger(ca_jet);
    
    if (tagged_jet == 0 ) continue;
    
    //Make sure jets have associated cluster sequence
    bool hasCS = tagged_jet.has_valid_cluster_sequence();
    if( !hasCS ) continue;

    higgs_candidates.push_back(tagged_jet);
  }

  //Require at least two mass-drop tagged large VR jets
  int const njet_cut=2;
  if((int)higgs_candidates.size() < njet_cut) 
  {
    Cut("Basic: Two MD FatJets",event_weight);
    event_weight=0;
    return;
  }

	//std::vector<double> split12_vec;
	//std::vector<double> tau21_vec;

	// Fails cuts
	if(event_weight<1e-30) return;
	
	//------------------------------------------------------------------
	// pt ordering is now done in the JetCluster_LargeVR(...) function
	//------------------------------------------------------------------

	// First of all, after basic selection, require that both jets are above 100 GeV
	double const pt_fatjet_ox = 100.0;
	// they should also be in central rapodity, |eta| < 2.5
	double const eta_fatjet_ox = 2.5;


	int const njet=2;
	// Restrict to the four leading jets in the event
	for(int ijet=0; ijet<njet;ijet++)
	{
		if(higgs_candidates.at(ijet).pt() < pt_fatjet_ox || 
			fabs( higgs_candidates.at(ijet).eta() ) > eta_fatjet_ox) 
			{
				Cut("Two dijets", event_weight);	// Kinematics cut on b-jets
				return;
			}
	}

	// Construct the Higgs candidates
	const fastjet::PseudoJet higgs1 = higgs_candidates.at(0);
	const fastjet::PseudoJet higgs2 = higgs_candidates.at(1);

	// Higgs mass window condition
	const double mass_diff1 = fabs(higgs1.m()-m_higgs)/m_higgs;
	const double mass_diff2 = fabs(higgs2.m()-m_higgs)/m_higgs;
	if( mass_diff1 > mass_resolution || mass_diff2 > mass_resolution ) 
	{
		Cut("Higgs window", event_weight);
		return;
	}
	// Histograms for the pt of the HH system
	// no cuts are applied on this variable
	const fastjet::PseudoJet dihiggs= higgs1+higgs2;


	// Now save the ntuples to be used by the TMVA or the ANNs
	//"# signal source m2fj pthh y2fj mHiggs1 mHiggs2 DeltaR_fj1fj2"
//	outputNTuple <<signal <<"\t"<<GetSample()<<"\t"<<dihiggs.m()<<"\t"<<dihiggs.pt()<<"\t"<<dihiggs.rapidity()<<"\t"<<
//	higgs1.m()<<"\t"<<higgs2.m()<<"\t"<<
//	split12_vec.at(0)<<"\t"<<split12_vec.at(1)<<"\t"<<
//	tau21_vec.at(0)<<"\t"<<tau21_vec.at(1)<<"\t"<<
//	fatjets.at(0).delta_R(fatjets.at(1))<<std::endl; 

	// Pass event
	Pass(event_weight);
}

/*
This routine read the event kinematics and performs the jet clustering
It also checkes that energy momentum is conserved event by event
This applies for small R jet clustering with the anti-kt algorithm
 */

void OxfordBoostVRAnalysis::JetCluster_LargeVR(fastjet::ClusterSequence const& cs_akt, std::vector<fastjet::PseudoJet>& higgs_candidates, double& event_weight)
{
  // Perform jet clustering with VR anti-kT
  // Note that here we use a small R clustering
 
  // Get all the jets (no pt cut here)
  std::vector<fastjet::PseudoJet> jets_vr_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
  VerifyFourMomentum(jets_vr_akt);
  
  // Calculate some substructure variables
  std::vector<double> split12_vec;
  std::vector<double> tau21_vec;
  split12_vec = SplittingScales( jets_vr_akt );
  tau21_vec = NSubjettiness( jets_vr_akt, jet_Rmax, jet_Rmin, jet_Rho );

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

  // By looking at the jet constituents
  // we can simulate the effects of b tagging
  const double initial_weight = event_weight;
  for(unsigned int ijet=0; ijet<jets_vr_akt.size(); ijet++)
  	if( TwoBTagging(jets_vr_akt[ijet]) )   // Check if at least two of its constituents are b quarks
  	{
  		higgs_candidates.push_back(jets_vr_akt.at(ijet));
  		event_weight *= btag_prob; // Account for b tagging efficiency
  	}
  	else // Else, account for the fake b-tag probabililty
  	{
  		higgs_candidates.push_back(jets_vr_akt.at(ijet));
  		event_weight *= btag_mistag;
  	}

  Cut("Basic: bTagging", initial_weight - event_weight);
  return;

} 


// ----------------------------------------------------------------------------------

bool OxfordBoostVRAnalysis::TwoBTagging( fastjet::PseudoJet const& jet ) const
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
	
	if(countB>1) return true;

 	return false; // no b-jets found
}
