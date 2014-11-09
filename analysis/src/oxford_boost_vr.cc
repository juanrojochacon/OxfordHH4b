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

OxfordBoostVRAnalysis::OxfordBoostVRAnalysis(std::string const& sampleName):
Analysis("oxford_boost_vr", sampleName)
{
	// Plotting parameters
	const double ptfj_min=0;
	const double ptfj_max=900;
	const int nbin_ptfj=20;

	const double DeltaRmin = 0;
	const double DeltaRmax = 5;
	
	// Fat Jet histograms (before kinematic cuts)
	BookHistogram(new YODA::Histo1D(nbin_ptfj, ptfj_min, ptfj_max), "ptfj1");
	BookHistogram(new YODA::Histo1D(nbin_ptfj, ptfj_min, ptfj_max), "ptfj2");

	BookHistogram(new YODA::Histo1D(20, 0, 200), "mfj1");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "mfj2");
	
	BookHistogram(new YODA::Histo1D(20, 0, 200), "split12_fj1");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "split12_fj2");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "tau21_fj1");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "tau21_fj2");
	
	// 2 fat jet system histograms
	BookHistogram(new YODA::Histo1D(20, 200, 1500), "m2fj");
	BookHistogram(new YODA::Histo1D(20, -2.5, 2.5), "y2fj");

	BookHistogram(new YODA::Histo1D(20, 0, 200), "mHiggs1");
	BookHistogram(new YODA::Histo1D(20, 0, 200), "mHiggs2");

	// Higgs histograms
	BookHistogram(new YODA::Histo1D(20, 0, 500), "pthh");
	BookHistogram(new YODA::Histo1D(20, 0, 600), "pth");

	BookHistogram(new YODA::Histo1D(20, DeltaRmin, DeltaRmax), "DeltaR_fj1fj2");

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
	std::vector<double> split12_vec;
	std::vector<double> tau21_vec;
	JetCluster_LargeVR(fs, fatjets, split12_vec, tau21_vec, event_weight);

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
		if(fatjets.at(ijet).pt() < pt_fatjet_ox || 
			fabs( fatjets.at(ijet).eta() ) > eta_fatjet_ox) 
			{
				Cut("Two dijets", event_weight);	// Kinematics cut on b-jets
				return;
			}
	}

	// Construct the Higgs candidates
	const fastjet::PseudoJet higgs1 = fatjets.at(0);
	const fastjet::PseudoJet higgs2 = fatjets.at(1);

	// Histograms before cut in pt Higgs candidates
	FillHistogram("pth", event_weight, higgs1.pt() );
	FillHistogram("pth", event_weight, higgs2.pt() );

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
	FillHistogram("pthh", event_weight, dihiggs.pt() );

	// Now save the ntuples to be used by the TMVA or the ANNs
	//"# signal source m2fj pthh y2fj mHiggs1 mHiggs2 DeltaR_fj1fj2"
	outputNTuple <<signal <<"\t"<<GetSample()<<"\t"<<dihiggs.m()<<"\t"<<dihiggs.pt()<<"\t"<<dihiggs.rapidity()<<"\t"<<
	higgs1.m()<<"\t"<<higgs2.m()<<"\t"<<
	split12_vec.at(0)<<"\t"<<split12_vec.at(1)<<"\t"<<
	tau21_vec.at(0)<<"\t"<<tau21_vec.at(1)<<"\t"<<
	fatjets.at(0).delta_R(fatjets.at(1))<<std::endl; 

	// Fill remaining histograms
	FillHistogram("m2fj", event_weight, dihiggs.m() );
	FillHistogram("y2fj", event_weight, dihiggs.rapidity() );

	FillHistogram("mHiggs1", event_weight, higgs1.m() );
	FillHistogram("mHiggs2", event_weight, higgs2.m() );

	FillHistogram("DeltaR_fj1fj2", event_weight, fatjets.at(0).delta_R(fatjets.at(1)) );

	// Pass event
	Pass(event_weight);
}

/*
This routine read the event kinematics and performs the jet clustering
It also checkes that energy momentum is conserved event by event
This applies for small R jet clustering with the anti-kt algorithm
 */

void OxfordBoostVRAnalysis::JetCluster_LargeVR(finalState const& particles, std::vector<fastjet::PseudoJet>& fatjets, std::vector<double>& split12_vec, std::vector<double>& tau21_vec, double& event_weight)
{
  // Perform jet clustering with VR anti-kT
  // Note that here we use a small R clustering
  
  // Set VR parameters
  static double const jet_Rmax	=1.0;
  static double const jet_Rmin	=0.1;
  static double const jet_Rho	=250.;
  
  // Set parameters for mass-drop tagger
  // mu = 0.67 and y = 0.09 are the default choice in FastJet
  double const mu = 0.67;
  double const ycut = 0.09;
  // The choice ycut = 0.15 is also recommended by Gavin to optimize S/sqrt(B)
  
  //Instantiate VR plugin
  VariableRPlugin lvjet_pluginAKT(jet_Rho, jet_Rmin, jet_Rmax, VariableRPlugin::AKTLIKE);
  fastjet::JetDefinition VR_AKT(&lvjet_pluginAKT);

  // Cluster all particles
  // The cluster sequence has to be saved to be used for jet substructure
  fastjet::ClusterSequence cs_akt(particles, VR_AKT);
  // Get all the jets (no pt cut here)
  std::vector<fastjet::PseudoJet> jets_vr_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
  VerifyFourMomentum(jets_vr_akt);
  
  // Calculate some substructure variables
  split12_vec = SplittingScales( jets_vr_akt );
  tau21_vec = NSubjettiness( jets_vr_akt, jet_Rmax, jet_Rmin, jet_Rho );

  // Fill the histograms for the pt of the fat jets before 
  // the corresponding kinematical cuts
  FillHistogram("ptfj1", event_weight, jets_vr_akt.at(0).pt() );
  FillHistogram("ptfj2", event_weight, jets_vr_akt.at(1).pt() );

  FillHistogram("mfj1", event_weight, jets_vr_akt.at(0).m() );
  FillHistogram("mfj2", event_weight, jets_vr_akt.at(1).m() );

  FillHistogram("split12_fj1", event_weight, split12_vec.at(0) );
  FillHistogram("split12_fj2", event_weight, split12_vec.at(1) );
  
  FillHistogram("tau21_fj1", event_weight, tau21_vec.at(0) );
  FillHistogram("tau21_fj2", event_weight, tau21_vec.at(1) );
  
  //-----------------------------------------------
  // Mass drop and b-tagging
  //-----------------------------------------------
  vector<fastjet::PseudoJet> higgs_candidates;

  unsigned int njet = jets_vr_akt.size();
  for (unsigned i = 0; i < njet; i++) {
  
    // first recluster with some large CA (needed for mass-drop)
    // Would it help to switch to CA VR to avoid this step?
    double const R_CA_boost = 1.2;
    fastjet::JetDefinition CA12(fastjet::cambridge_algorithm, R_CA_boost);
    fastjet::ClusterSequence cs_sub(jets_vr_akt[i].constituents(), CA12);
    
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

  // By looking at the jet constituents
  // we can simulate the effects of b tagging
  
  // Here we would need to implement double b-tagging to 
  // improve background rejection

  // Loop over the 2 hardest jets in event only
  const double initial_weight = event_weight;
  unsigned int nhjet = higgs_candidates.size();
  
  for(unsigned int ijet=0; ijet<nhjet; ijet++)
  {
	if( TwoBTagging(higgs_candidates[ijet]) )   // Check if at least two of its constituents are b quarks
	{
		fatjets.push_back(higgs_candidates.at(ijet));
		event_weight *= btag_prob; // Account for b tagging efficiency
	}
	else // Else, account for the fake b-tag probabililty
	{
		fatjets.push_back(higgs_candidates.at(ijet));
		event_weight *= btag_mistag;
	}

	// cut from btagging
	Cut("Basic: bTagging", initial_weight - event_weight);
  }
} 

// ----------------------------------------------------------------------------------

bool OxfordBoostVRAnalysis::BTagging( fastjet::PseudoJet const& jet ) const
{
	// Cuts for the b-jet candidates for b-tagging
	double const pt_btagging=15;
	
	bool hasCS = jet.has_valid_cluster_sequence();
	if( !hasCS ) return false;

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

// ----------------------------------------------------------------------------------
// Recluster with kt algorithm to obtain splitting scales
std::vector< double > OxfordBoostVRAnalysis::SplittingScales( std::vector<fastjet::PseudoJet> jetVec )
{
   
   //vectors that contain the respective splitting scales for all jets
   std::vector<double> split12_vec;
   
   for( int i = 0; i < (int) jetVec.size(); i++){
   
      // For now: Calculate substructure information only for the two leading jets
      if( i > 1 ) continue;
   
      double split12 = -1.;

      if (!jetVec.at(i).has_constituents()){
         std::cout << "ERROR! Splittings (d12, d23, ...) can only be calculated on jets for which the constituents are known."<< std::endl;
         split12_vec.push_back(-1);
         continue;
      }
   
      vector<fastjet::PseudoJet> constits = jetVec.at(i).constituents();
      
//       std::cout << "Jet " <<  i << " has " << constits.size() << " constituents." << std::endl;
      
      fastjet::JetDefinition ekt_jd = fastjet::JetDefinition( fastjet::kt_algorithm, 1.5, fastjet::E_scheme, fastjet::Best);
      const fastjet::ClusterSequence kt_seq_excl = fastjet::ClusterSequence( constits, ekt_jd);
      fastjet::PseudoJet kt_jet = sorted_by_pt( kt_seq_excl.inclusive_jets())[0];
      
      split12 = 1.5*sqrt( kt_seq_excl.exclusive_subdmerge( kt_jet, 1));
      
      split12_vec.push_back(split12);
   }
   
   return split12_vec;
}

// ----------------------------------------------------------------------------------
// Recluster with kt algorithm to obtain nsubjettiness
std::vector< double > OxfordBoostVRAnalysis::NSubjettiness( std::vector<fastjet::PseudoJet> jetVec, double jet_Rmax, double jet_Rmin, double jet_Rho )
{

   // Reclustering with VR kt algorithm to obtain nsubjettiness
   
   //vector that contain the respective nsubjettiness variables for all jets
   std::vector<double> tau1_vec;
   std::vector<double> tau2_vec;
   std::vector<double> tau3_vec;
   
   std::vector<double> tau21_vec;
   
   double alpha=1;
   
   for( int i = 0; i < (int) jetVec.size(); i++){
   
      // For now: Calculate substructure information only for the two leading jets
      if( i > 1 ) continue;
   
      double tau1 = -1.;
      double tau2 = -1.;
      double tau3 = -1.;
      
      // Calculate effective jet radius
      double jet_Pt = jetVec[i].pt();
      
      double jet_rad;
      if( jet_Pt > jet_Rmax ) jet_rad = jet_Rmax;
      else if( jet_Pt < jet_Rmin ) jet_rad = jet_Rmin;
      else jet_rad = jet_Rho / jet_Pt;

      if (!jetVec.at(i).has_constituents()){
         std::cout << "ERROR! NSubjettiness (tau1, tau2, ...) can only be calculated on jets for which the constituents are known."<< std::endl;
         
         tau1_vec.push_back(-1);
	 tau2_vec.push_back(-1);
	 tau3_vec.push_back(-1);
	
	 double tau21 = tau2/tau1;
	 tau21_vec.push_back(-1);
         continue;
      }
      
      //Code snippet taken from https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetSubstructureVariables#N_subjettiness
      vector<fastjet::PseudoJet> constits = jetVec.at(i).constituents();
      if(constits.size()>0){
      
	 double R_kt = 1.5;
	 fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::kt_algorithm,R_kt,fastjet::E_scheme,fastjet::Best);
	 fastjet::ClusterSequence kt_clust_seq(constits, jet_def);
	 vector<fastjet::PseudoJet> kt1axes = kt_clust_seq.exclusive_jets(1);
	 double tauNum = 0.0;
	 double tauDen = 0.0;
	 for (int i = 0; i < (int)constits.size(); i++) {
	    // find minimum distance
	    double minR = 10000.0; // large number
	    for (int j = 0; j < (int)kt1axes.size(); j++) {
	       double tempR = sqrt(constits[i].squared_distance(kt1axes[j])); // delta R distance
	       if (tempR < minR) minR = tempR;
	    }
	    tauNum += constits[i].perp() * pow(minR,alpha);
	    tauDen += constits[i].perp() * pow(jet_rad,alpha);
	 }
	 tau1 = tauNum/tauDen;

	 if(constits.size()>1){
	    vector<fastjet::PseudoJet> kt2axes = kt_clust_seq.exclusive_jets(2);
	    tauNum = 0.0;
	    tauDen = 0.0;
	    for (int i = 0; i < (int)constits.size(); i++) {
	    // find minimum distance
	    double minR = 10000.0; // large number
	    for (int j = 0; j < (int)kt2axes.size(); j++) {
	       double tempR = sqrt(constits[i].squared_distance(kt2axes[j]));
	       if (tempR < minR) minR = tempR;
	    }
	    tauNum += constits[i].perp() * pow(minR,alpha);
	    tauDen += constits[i].perp() * pow(jet_rad,alpha);
	    }
	    tau2 = tauNum/tauDen;
	    
	    if(constits.size() > 2){
	    
	       vector<fastjet::PseudoJet> kt3axes = kt_clust_seq.exclusive_jets(3);
	       tauNum = 0.0;
	       tauDen = 0.0;
	       for (int i = 0; i < (int)constits.size(); i++) {
	       // find minimum distance
	       double minR = 10000.0; // large number
	       for (int j = 0; j < (int)kt3axes.size(); j++) {
		  double tempR = sqrt(constits[i].squared_distance(kt3axes[j]));
		  if (tempR < minR) minR = tempR;
	       }
	       tauNum += constits[i].perp() * pow(minR,alpha);
	       tauDen += constits[i].perp() * pow(jet_rad,alpha);
	       }
	       tau3 = tauNum/tauDen;
	    }

	 }
      }
      
      tau1_vec.push_back(tau1);
      tau2_vec.push_back(tau2);
      tau3_vec.push_back(tau3);
      
      double tau21 = tau2/tau1;
      tau21_vec.push_back(tau21);
   }
   
   return tau21_vec;
}
