// oxford_boost_fr.cc

#include "oxford_boost_fr.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

#include "YODA/Histo1D.h"

using namespace fastjet::contrib;

// Perform jet clustering with anti-kT
// Note that here we use a small R clustering
static double const jetR=1.0; // To avoid overlapping b's as much as possible
fastjet::JetDefinition akt(fastjet::antikt_algorithm, jetR);

// first recluster again with the large-R but with Cambridge-Aachen
static double const jetR_1p2=1.2; // To try to merge two b quarks into the same jet
static fastjet::JetDefinition CA10(fastjet::cambridge_algorithm, jetR_1p2);
// Set parameters of the mass drop tagger for jet substructure
// mu = 0.67 and y = 0.09 are the default choice in FastJet
static const double mu = 0.67;
static const double ycut = 0.09;


OxfordBoostFRAnalysis::OxfordBoostFRAnalysis(std::string const& sampleName):
Analysis("oxford_boost_fr", sampleName)
{
  // Plotting parameters
  const double ptfj_min=0;
  const double ptfj_max=900;

  const double DeltaRmin = 1;
  const double DeltaRmax = 5;
  
  // *********************** preCut **************************

  // Fat Jet histograms
  BookHistogram(new YODA::Histo1D(20, 0, 20), "nFatJets_preCut");
  
  BookHistogram(new YODA::Histo1D(30, ptfj_min, ptfj_max), "ptfj1_preCut");
  BookHistogram(new YODA::Histo1D(30, ptfj_min, ptfj_max), "ptfj2_preCut");

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
  
  // ************************* before b-tagging ********************************
  // Subjet histograms (after pt, m cuts on both fatjets)
  BookHistogram(new YODA::Histo1D(20, 0, 20), "nSubjets_fj0_preBtag");
  BookHistogram(new YODA::Histo1D(20, 0, 20), "nSubjets_fj1_preBtag");

  BookHistogram(new YODA::Histo1D(20, 0, 500), "leadSubjet_fj0_pt_preBtag");
  BookHistogram(new YODA::Histo1D(20, 0, 500), "leadSubjet_fj1_pt_preBtag");
  
  BookHistogram(new YODA::Histo1D(20, 0, 500), "subleadSubjet_fj0_pt_preBtag");
  BookHistogram(new YODA::Histo1D(20, 0, 500), "subleadSubjet_fj1_pt_preBtag");

  BookHistogram(new YODA::Histo1D(10, 0, 10), "leadSubjet_fj0_nBQuarks_preBtag");
  BookHistogram(new YODA::Histo1D(10, 0, 10), "leadSubjet_fj1_nBQuarks_preBtag");

  BookHistogram(new YODA::Histo1D(10, 0, 10), "subleadSubjet_fj0_nBQuarks_preBtag");
  BookHistogram(new YODA::Histo1D(10, 0, 10), "subleadSubjet_fj1_nBQuarks_preBtag");
  
  BookHistogram(new YODA::Histo1D(10, 0, 10), "leadSubjet_fj0_nCQuarks_preBtag");
  BookHistogram(new YODA::Histo1D(10, 0, 10), "leadSubjet_fj1_nCQuarks_preBtag");
  
  BookHistogram(new YODA::Histo1D(10, 0, 10), "subleadSubjet_fj0_nCQuarks_preBtag");
  BookHistogram(new YODA::Histo1D(10, 0, 10), "subleadSubjet_fj1_nCQuarks_preBtag");
  
  
  // ************************* after b-tagging ********************************
  BookHistogram(new YODA::Histo1D(20, 0, 20), "nSubjets_fj0_postBtag");
  BookHistogram(new YODA::Histo1D(20, 0, 20), "nSubjets_fj1_postBtag");

  BookHistogram(new YODA::Histo1D(20, 0, 500), "leadSubjet_fj0_pt_postBtag");
  BookHistogram(new YODA::Histo1D(20, 0, 500), "leadSubjet_fj1_pt_postBtag");
  
  BookHistogram(new YODA::Histo1D(20, 0, 500), "subleadSubjet_fj0_pt_postBtag");
  BookHistogram(new YODA::Histo1D(20, 0, 500), "subleadSubjet_fj1_pt_postBtag");

  BookHistogram(new YODA::Histo1D(10, 0, 10), "leadSubjet_fj0_nBQuarks_postBtag");
  BookHistogram(new YODA::Histo1D(10, 0, 10), "leadSubjet_fj1_nBQuarks_postBtag");

  BookHistogram(new YODA::Histo1D(10, 0, 10), "subleadSubjet_fj0_nBQuarks_postBtag");
  BookHistogram(new YODA::Histo1D(10, 0, 10), "subleadSubjet_fj1_nBQuarks_postBtag");
  
  BookHistogram(new YODA::Histo1D(10, 0, 10), "leadSubjet_fj0_nCQuarks_postBtag");
  BookHistogram(new YODA::Histo1D(10, 0, 10), "leadSubjet_fj1_nCQuarks_postBtag");
  
  BookHistogram(new YODA::Histo1D(10, 0, 10), "subleadSubjet_fj0_nCQuarks_postBtag");
  BookHistogram(new YODA::Histo1D(10, 0, 10), "subleadSubjet_fj1_nCQuarks_postBtag");
  
  // ************************* postCut ********************************

  // Fat Jet histograms
  BookHistogram(new YODA::Histo1D(20, 0, 20), "nFatJets_postCut");
  
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

	const std::string tupleSpec = "# signal source weight m2fj pthh y2fj mHiggs1 mHiggs2 split12_Higgs1 split12_Higgs2 tau21_Higgs1 tau21_Higgs2 DeltaR_fj1fj2";
	outputNTuple<<tupleSpec<<std::endl;

  // Order cutflow
  Cut("Basic: Fatjet kinematic cuts ", 0);
  Cut("Basic: 2 subjets for each fatjet ",0);
  Cut("Basic: bTagging", 0);
  Cut("Basic: No double b-tagged subjets ",0);
  Cut("BDRS mass-drop", 0);
  //Cut("fatjet deltaEta", 0);
  Cut("Higgs window", 0);
}

void OxfordBoostFRAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
  Analysis::Analyse(signal, weightnorm, fs);

	// Set initial weight
	double event_weight = weightnorm;

	// Perform jet clustering with VR anti-kT
	std::vector<fastjet::PseudoJet> fatjets;
	
	// Fetch jets and substructure information
	std::vector<double> split12_vec;
	std::vector<double> tau21_vec;
	JetCluster_LargeFR(fs, fatjets, split12_vec, tau21_vec, event_weight);

	// Fails cuts
	if(event_weight<1e-30) return;// Cut("Rounding", event_weight);
	
	//------------------------------------------------------------------------------
	// pt ordering and basic cuts now done in the JetCluster_LargeFR(...) function
	//------------------------------------------------------------------------------

	// Construct the Higgs candidates
	const fastjet::PseudoJet higgs1 = fatjets.at(0);
	const fastjet::PseudoJet higgs2 = fatjets.at(1);

	// Histograms for the HH system
	const fastjet::PseudoJet dihiggs= higgs1+higgs2;


// ************************************************************************************

  // Same as in the UCL analysis
  // require that these two leading fat jets are not too separated in rapidity
  //const double delta_eta_dijet_fatjet=1.5;
  //const double delta_eta_dijet = fabs(fatjets.at(0).eta()- fatjets.at(1).eta());
  //if(delta_eta_dijet > delta_eta_dijet_fatjet) return Cut("fatjet deltaEta", event_weight);

  // Higgs mass window condition
  const double mass_diff1 = fabs(fatjets.at(0).m()-m_higgs)/m_higgs;
  const double mass_diff2 = fabs(fatjets.at(1).m()-m_higgs)/m_higgs;
  if( mass_diff1 > mass_resolution || mass_diff2 > mass_resolution ) 
    return  Cut("Higgs window", event_weight);


// ************************************************************************************


	int nJets = fatjets.size();
	// Record jet multiplicity after cuts
	FillHistogram("nFatJets_postCut", event_weight, nJets );

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
	outputNTuple <<signal <<"\t"<<GetSample()<<"\t"<<event_weight<<"\t"<<dihiggs.m()<<"\t"<<dihiggs.pt()<<"\t"<<dihiggs.rapidity()<<"\t"<<
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

void OxfordBoostFRAnalysis::JetCluster_LargeFR(finalState const& fs, std::vector<fastjet::PseudoJet>& fatjets, std::vector<double>& split12_vec, std::vector<double>& tau21_vec, double& event_weight)
{
  
  finalState fsc;
  //Select only charged fs particles
  for(int i=0; i<(int)fs.size(); i++){
          int userid = fs.at(i).user_index();
          //std::cout << "userid " << userid << std::endl;
          if(abs(userid)<6) fsc.push_back(fs.at(i));
  }

  // Get all the jets (no pt cut here)
  fastjet::ClusterSequence cs_akt(fs, akt);
  std::vector<fastjet::PseudoJet> jets_fr_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
  //VerifyFourMomentum(jets_fr_akt);
  
  int nJets = jets_fr_akt.size();
  // We require at least 2 fatjets in the event, else discard event
  int const njet=2;
  if(nJets < njet) 
  {
	  Cut("Basic: Two fatjets",event_weight);
	  event_weight=0;
	  return;
  }
  
  // Calculate some substructure variables
  split12_vec = SplittingScales( jets_fr_akt );
  tau21_vec = NSubjettiness( jets_fr_akt, jetR );
  
  const fastjet::PseudoJet dihiggs= jets_fr_akt.at(0)+jets_fr_akt.at(1);
  
  // Record jet multiplicity before cuts
  FillHistogram("nFatJets_preCut", event_weight, nJets );

  // Fill the histograms for the pt of the fat jets before 
  // the corresponding kinematical cuts
  FillHistogram("ptfj1_preCut", event_weight, jets_fr_akt.at(0).pt() );
  FillHistogram("ptfj2_preCut", event_weight, jets_fr_akt.at(1).pt() );

  FillHistogram("mfj1_preCut", event_weight, jets_fr_akt.at(0).m() );
  FillHistogram("mfj2_preCut", event_weight, jets_fr_akt.at(1).m() );

  FillHistogram("split12_fj1_preCut", event_weight, split12_vec.at(0) );
  FillHistogram("split12_fj2_preCut", event_weight, split12_vec.at(1) );
  
  FillHistogram("tau21_fj1_preCut", event_weight, tau21_vec.at(0) );
  FillHistogram("tau21_fj2_preCut", event_weight, tau21_vec.at(1) );
  
  FillHistogram("m2fj_preCut", event_weight, dihiggs.m() );
  FillHistogram("pT2fj_preCut", event_weight, dihiggs.pt() );
  FillHistogram("y2fj_preCut", event_weight, dihiggs.rapidity() );
  FillHistogram("DeltaR_fj1fj2_preCut", event_weight, jets_fr_akt.at(0).delta_R(jets_fr_akt.at(1)) );
  
  
  //------------------------------------------------------------------------------------------------------
  // First of all, after basic selection, require that both jets are above 100 GeV
  double const pt_fatjet_ox = 250.0;
  // they should also be in central rapidity, |eta| < 2.5
  double const eta_fatjet_ox = 2.5;
  // they should have a mass in a window around m_H
  double const m_fatjet_ox_upper = 160.;
  double const m_fatjet_ox_lower = 90.;
  
  //Basic kinematic cuts on the jets
  int const nfatjet=2;
  // Restrict to the two leading jets in the event
  for(int ijet=0; ijet<nfatjet;ijet++)
  {
	  if(jets_fr_akt.at(ijet).pt() < pt_fatjet_ox || 
		  fabs( jets_fr_akt.at(ijet).eta() ) > eta_fatjet_ox)
		  {
		      // Higgs mass window cut
		      if( jets_fr_akt.at(ijet).m() > m_fatjet_ox_upper || jets_fr_akt.at(ijet).m() < m_fatjet_ox_lower )
		      {
	                  Cut("Basic: Fatjet kinematic cuts ", event_weight);
			  event_weight=0;
			  return;
		      }
		  }
		  else fatjets.push_back(jets_fr_akt.at(ijet));
  }
  
  //------------------------------------------------------------------------------------------------------
  // Simulate subjets:
  // Perform small-R jet clustering with anti-kT
  static double const jetR=0.3; // To avoid overlapping b's as much as possible
  fastjet::JetDefinition jd_subjets(fastjet::antikt_algorithm, jetR);
  fastjet::ClusterSequence cs_subjets(fsc, jd_subjets);

  std::vector<fastjet::PseudoJet> jets_akt = sorted_by_pt( cs_subjets.inclusive_jets()  );
  
  std::vector<fastjet::PseudoJet> subjets;

  //---------------------------------------------------------------------------------------
  // Ghost association of small jets to fat jets
  std::vector<fastjet::PseudoJet> subjets_fj0_unsrt;
  std::vector<fastjet::PseudoJet> subjets_fj1_unsrt;
  
  // This function considers only small-R jets with pT > 20 GeV and |eta| < 2.5 for ghost association
  get_assoc_trkjets( jets_fr_akt.at(0), jets_akt, subjets_fj0_unsrt, false);
  get_assoc_trkjets( jets_fr_akt.at(1), jets_akt, subjets_fj1_unsrt, false);
  
  // Subjet multiplicity (before any subjet cuts)
  FillHistogram("nSubjets_fj0_preBtag", event_weight, (int)subjets_fj0_unsrt.size() );
  FillHistogram("nSubjets_fj1_preBtag", event_weight, (int)subjets_fj1_unsrt.size() );
  
  // Require at least 2 subjets per fat jet
  int const nsubjet=2;
  if( ((int)subjets_fj0_unsrt.size() < nsubjet ) || ((int)subjets_fj1_unsrt.size() < nsubjet ) ) 
  {
          Cut("Basic: 2 subjets for each fatjet ", event_weight);
	  event_weight=0;
	  return;
  }
  
  std::vector<fastjet::PseudoJet> subjets_fj0 = sorted_by_pt( subjets_fj0_unsrt  );
  std::vector<fastjet::PseudoJet> subjets_fj1 = sorted_by_pt( subjets_fj1_unsrt  );
  
  FillHistogram("leadSubjet_fj0_pt_preBtag", event_weight, subjets_fj0[0].pt() );
  FillHistogram("leadSubjet_fj1_pt_preBtag", event_weight, subjets_fj1[0].pt() );
  FillHistogram("subleadSubjet_fj0_pt_preBtag", event_weight, subjets_fj0[1].pt() );
  FillHistogram("subleadSubjet_fj1_pt_preBtag", event_weight, subjets_fj1[1].pt() );

  int nBQuarks_fj0_0 = BTagging(subjets_fj0[0]);
  int nBQuarks_fj0_1 = BTagging(subjets_fj0[1]);
  int nBQuarks_fj1_0 = BTagging(subjets_fj1[0]);
  int nBQuarks_fj1_1 = BTagging(subjets_fj1[1]);

  FillHistogram("leadSubjet_fj0_nBQuarks_preBtag", event_weight, nBQuarks_fj0_0 );
  FillHistogram("subleadSubjet_fj0_nBQuarks_preBtag", event_weight, nBQuarks_fj0_1 );
  FillHistogram("leadSubjet_fj1_nBQuarks_preBtag", event_weight, nBQuarks_fj1_0 );
  FillHistogram("subleadSubjet_fj1_nBQuarks_preBtag", event_weight, nBQuarks_fj1_1 );


  int nCQuarks_fj0_0 = CTagging(subjets_fj0[0]);
  int nCQuarks_fj0_1 = CTagging(subjets_fj0[1]);
  int nCQuarks_fj1_0 = CTagging(subjets_fj1[0]);
  int nCQuarks_fj1_1 = CTagging(subjets_fj1[1]);
  
  FillHistogram("leadSubjet_fj0_nCQuarks_preBtag", event_weight, nCQuarks_fj0_0 );
  FillHistogram("subleadSubjet_fj0_nCQuarks_preBtag", event_weight, nCQuarks_fj0_1 );
  FillHistogram("leadSubjet_fj1_nCQuarks_preBtag", event_weight, nCQuarks_fj1_0 );
  FillHistogram("subleadSubjet_fj1_nCQuarks_preBtag", event_weight, nCQuarks_fj1_1 );
  //---------------------------------------------------------------------------------------
  
  std::vector<fastjet::PseudoJet> bjets_jet0;
  std::vector<fastjet::PseudoJet> bjets_jet1;
  
  // By looking at the jet constituents
  // we can simulate the effects of b tagging
  // Here I am only attempting to b-tag the hardest two jets
  const double initial_weight = event_weight;
  for(unsigned int ijet=0; ijet<2; ijet++)
  {
    bjets_jet0.push_back(subjets_fj0.at(ijet));

  	if( BTagging(subjets_fj0[ijet]) )   // Check if at least two of its constituents are b quarks
  	{
		const double btag_prob = btag_eff( subjets_fj0.at(ijet).pt() );
  		event_weight *= btag_prob; // Account for b tagging efficiency
  	}
  	else if( CTagging(subjets_fj0[ijet]) ) 
	{
  		const double ctag_prob = charm_eff( subjets_fj0.at(ijet).pt() );
  		event_weight *= ctag_prob; // Account for c (mis-)tagging efficiency
  	}
  	else // Else, account for the fake b-tag probability
  	{
		const double mistag_prob = mistag_eff( subjets_fj0.at(ijet).pt() );
  		event_weight *= mistag_prob;
  	}
  }

  for(unsigned int ijet=0; ijet<2; ijet++)
  {
    bjets_jet1.push_back(subjets_fj1.at(ijet));

  	if( BTagging(subjets_fj1[ijet]) )   // Check if at least one of its constituents are b quarks
  	{
  		const double btag_prob = btag_eff( subjets_fj1.at(ijet).pt() );
  		event_weight *= btag_prob; // Account for b tagging efficiency
  	}
  	else if( CTagging(subjets_fj1[ijet]) ) 
    {
  		const double ctag_prob = charm_eff( subjets_fj1.at(ijet).pt() );
  		event_weight *= ctag_prob; // Account for c (mis-)tagging efficiency
  	}
  	else // Else, account for the fake b-tag probability
  	{	
		  const double mistag_prob = mistag_eff( subjets_fj1.at(ijet).pt() );
  		event_weight *= mistag_prob;
  	}
  }
  
  Cut("Basic: bTagging", initial_weight - event_weight);
  
  
  //---------------------------------------------------------------------------------------
  // Fill post-b-tagging histograms
  FillHistogram("leadSubjet_fj0_pt_postBtag", event_weight, subjets_fj0[0].pt() );
  FillHistogram("leadSubjet_fj1_pt_postBtag", event_weight, subjets_fj1[0].pt() );
  FillHistogram("subleadSubjet_fj0_pt_postBtag", event_weight, subjets_fj0[1].pt() );
  FillHistogram("subleadSubjet_fj1_pt_postBtag", event_weight, subjets_fj1[1].pt() );
  
  FillHistogram("leadSubjet_fj0_nBQuarks_postBtag", event_weight, nBQuarks_fj0_0 );
  FillHistogram("subleadSubjet_fj0_nBQuarks_postBtag", event_weight, nBQuarks_fj0_1 );
  FillHistogram("leadSubjet_fj1_nBQuarks_postBtag", event_weight, nBQuarks_fj1_0 );
  FillHistogram("subleadSubjet_fj1_nBQuarks_postBtag", event_weight, nBQuarks_fj1_1 );
  
  FillHistogram("leadSubjet_fj0_nCQuarks_postBtag", event_weight, nCQuarks_fj0_0 );
  FillHistogram("subleadSubjet_fj0_nCQuarks_postBtag", event_weight, nCQuarks_fj0_1 );
  FillHistogram("leadSubjet_fj1_nCQuarks_postBtag", event_weight, nCQuarks_fj1_0 );
  FillHistogram("subleadSubjet_fj1_nCQuarks_postBtag", event_weight, nCQuarks_fj1_1 );
  //---------------------------------------------------------------------------------------
  // TESTING DOUBLE BTagging
  int nbbJets = 0;
  
  for(int ijet=0; ijet<2; ijet++)
  {
	  int bQuarks = BTagging(subjets_fj0.at(ijet));

	  const double dice = ((double) rand() / (double)(RAND_MAX));
	  if( bQuarks > 1 )   // Check if at least two of its constituents are b quarks
	  {
		  if (dice < bbtag_prob)
			  nbbJets++;				
	  }
	  else if( bQuarks == 1) // Else, account for the fake bb-tag probabililty
	  {
		  if (dice < bbtag_mistag)
			  nbbJets++;
	  }
  }

  for(int ijet=0; ijet<2; ijet++)
  {
	  int bQuarks = BTagging(subjets_fj1.at(ijet));

	  const double dice = ((double) rand() / (double)(RAND_MAX));
	  if( bQuarks > 1 )   // Check if at least two of its constituents are b quarks
	  {
		  if (dice < bbtag_prob)
			  nbbJets++;				
	  }
	  else if( bQuarks == 1) // Else, account for the fake bb-tag probabililty
	  {
		  if (dice < bbtag_mistag)
			  nbbJets++;
	  }
  }
  
  // Return if there is a bb-tagged jet
  if(nbbJets >  0)
  {
	   Cut("Basic: No double b-tagged subjets ", event_weight);
     event_weight = 0;
     return;
  }
  

    // Now look for substructure in each of these two dijets using the BDRS mass-drop tagger
  int nTagged = 0;
  for (int i = 0; i < 2; i++) 
  {
    fastjet::ClusterSequence cs_sub(fatjets[i].constituents(), CA10);
    fastjet::PseudoJet ca_jet = sorted_by_pt(cs_sub.inclusive_jets())[0];

    // now run mass drop tagger
    // parameters are specified in settings.h
    fastjet::MassDropTagger md_tagger(mu, ycut);
    fastjet::PseudoJet tagged_jet = md_tagger(ca_jet);

    // If tagging succesful - declare as Higgs candidate
    if (tagged_jet != 0 )  
      nTagged++; 
  }

  // If we don't have a mass-drop tag in each of the two leading large-R jets
  // discard the event
  if(nTagged!=2) 
  {
    Cut("BDRS mass-drop", event_weight);
    event_weight = 0;
    return;
  }
  
  return;
}


// ----------------------------------------------------------------------------------

int OxfordBoostFRAnalysis::BTagging( fastjet::PseudoJet const& jet ) const
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


int OxfordBoostFRAnalysis::CTagging( fastjet::PseudoJet const& jet ) const
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
