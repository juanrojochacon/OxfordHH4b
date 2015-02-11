// oxford_res_fr.cc

#include "btag_test.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

#include "YODA/Histo1D.h"

using namespace fastjet::contrib;

// Jet definition
static double const jetR=0.4; // To avoid overlapping b's as much as possible
static const fastjet::JetDefinition akt(fastjet::antikt_algorithm, jetR);

// b tagging
// Choose working point with high purity
static double const test_btag_prob = 0.80; 		// Probability of correct b tagging
static double const test_btag_mistag = 0.01; 	// Mistag probability  
static double const test_bbtag_prob = 0.80; 		// Probability of correct b tagging
static double const test_bbtag_mistag = 0.01; 	// Mistag probability  

bTagTestAnalysis::bTagTestAnalysis(std::string const& sampleName):
Analysis("bTagTest_hardest4", sampleName)
{
	const std::string tupleSpec = "# signal source weight";
	outputNTuple<<tupleSpec<<std::endl;

	BookHistogram(new YODA::Histo1D(5, 0, 5), "truth_NbJets");
	BookHistogram(new YODA::Histo1D(5, 0, 5), "truth_NbConstituents");

	// Order cutflow
	Cut("Basic: Two dijets", 0);
	Cut("Basic: bTagging", 0);
	Cut("Basic: bbTagging", 0);
}

void bTagTestAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
	Analysis::Analyse(signal, weightnorm, fs);

	// Set initial weight
	double event_weight = weightnorm;

	// Fetch jets
	fastjet::ClusterSequence cs_akt(fs, akt);
	std::vector<fastjet::PseudoJet> jets_fr_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
	VerifyFourMomentum(jets_fr_akt);

	// We require at least 4 jets in the event, else discard event
	int const njet=4;
	if((int)jets_fr_akt.size() < njet) 
		return Cut("Basic: Two dijets",event_weight);

	// Determine b-jets
	const double initial_weight = event_weight;
	std::vector<fastjet::PseudoJet> bjets;

	// Loop over the 4 hardest jets in event only
	int NbJets = 0;
	for(int ijet=0; ijet<njet;ijet++)
	{
		int bQuarks = BTagging(jets_fr_akt[ijet]);
		FillHistogram("truth_NbConstituents", 1, bQuarks+0.5 );

		// Add jet
		bjets.push_back(jets_fr_akt.at(ijet));

		if( bQuarks > 0 )   // Check if at least one of its constituents is a b quark
			NbJets++;
	}

	FillHistogram("truth_NbJets", 1, NbJets+0.5 );

	// cut from btagging
	event_weight*=pow(test_btag_mistag, (double) (4-NbJets))*pow(test_btag_prob, (double)NbJets);
	Cut("Basic: bTagging", initial_weight - event_weight);
	
	
	// TESTING DOUBLE BTagging
	int nbbJets = 0;
	
	for(int ijet=0; ijet<njet; ijet++)
	{
		int bQuarks = BTagging(jets_fr_akt[ijet]);

		const double dice = ((double) rand() / (double)(RAND_MAX));
		if( bQuarks > 1 )   // Check if at least two of its constituents are b quarks
		{
			if (dice < test_bbtag_prob)
				nbbJets++;				
		}
		else if( bQuarks == 1) // Else, account for the fake bb-tag probabililty
		{
			if (dice < test_bbtag_mistag)
				nbbJets++;
		}
	}

	// Return if there is a bb-tagged jet
	if(nbbJets >  0)
		return 	Cut("Basic: bbTagging", event_weight);

	
	// ************************************* MVA Output **********************************************************
	outputNTuple <<signal <<"\t"<< GetSample() <<"\t"<<event_weight<<"\t"<<std::endl; 

	// Pass event
	Pass(event_weight);

}

// ----------------------------------------------------------------------------------

int bTagTestAnalysis::BTagging( fastjet::PseudoJet const& jet ) const
{
	// Cuts for the b-jet candidates for b-tagging
	double const pt_btagging=15;

	// Get the jet constituents
	const std::vector<fastjet::PseudoJet>& jet_constituents = jet.constituents();

	// Loop over constituents and look for b quarks
	// also b quarks must be above some minimum pt
	int bquarks = 0;
	for(size_t i=0; i<jet_constituents.size(); i++)
	{
		// Flavour of jet constituent
		const int userid= jet_constituents.at(i).user_index();
		const double pt_bcandidate = jet_constituents.at(i).pt();

		if(abs(userid) ==5 )
			if( pt_bcandidate > pt_btagging)
		  		bquarks++;
	}

 	return bquarks; // no b-jets found
}

// ****************************************** UCL Style b-Tagging test *********************************************

bTagTestUCLAnalysis::bTagTestUCLAnalysis(std::string const& sampleName):
Analysis("bTagTest_UCL", sampleName)
{
	const std::string tupleSpec = "# signal source weight";
	outputNTuple<<tupleSpec<<std::endl;

	BookHistogram(new YODA::Histo1D(5, 0, 5), "UCL_truth_NbJets");
	BookHistogram(new YODA::Histo1D(5, 0, 5), "UCL_truth_NbConstituents");

	// Order cutflow
	Cut("Basic: Two dijets pT > 40 GeV", 0);
	Cut("Basic: bTagging", 0);
}

void bTagTestUCLAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
	Analysis::Analyse(signal, weightnorm, fs);

	// Set initial weight
	double event_weight = weightnorm;

	// Fetch jets
	fastjet::ClusterSequence cs_akt(fs, akt);
	std::vector<fastjet::PseudoJet> jets_fr_akt = sorted_by_pt( cs_akt.inclusive_jets( 40.0 )  );
	//VerifyFourMomentum(jets_fr_akt); // not a good idea here!

	// We require at least 4 jets in the event, else discard event
	int const njet=4;
	if((int)jets_fr_akt.size() < njet) 
		return Cut("Basic: Two dijets pT > 40 GeV",event_weight);

	// Determine b-jets
	std::vector<fastjet::PseudoJet> bjets;

	// Loop over the 4 hardest jets in event only
	int NbJets = 0;
	for(int ijet=0; ijet<njet;ijet++)
	{
		int bQuarks = BTagging(jets_fr_akt[ijet]);
		FillHistogram("UCL_truth_NbConstituents", 1, bQuarks+0.5 );

		const double dice = ((double) rand() / (double)(RAND_MAX));
		if( bQuarks > 0 )   // Check if at least one of its constituents is a b quark
		{
			NbJets++;

			if (dice < test_btag_prob)
				bjets.push_back(jets_fr_akt.at(ijet));				
		}
		else // Else, account for the fake b-tag probabililty
		{
			if (dice < test_btag_mistag)
				bjets.push_back(jets_fr_akt.at(ijet));
		}
	}

	// Not enough b-jets
	if (bjets.size() < 4)
		return 	Cut("Basic: bTagging", event_weight);

	FillHistogram("UCL_truth_NbJets", 1, NbJets+0.5 );
	
	// ************************************* MVA Output **********************************************************
	outputNTuple <<signal <<"\t"<< GetSample() <<"\t"<<event_weight<<"\t"<<std::endl; 

	// Pass event
	Pass(event_weight);

}

// ----------------------------------------------------------------------------------

int bTagTestUCLAnalysis::BTagging( fastjet::PseudoJet const& jet ) const
{
	// Cuts for the b-jet candidates for b-tagging
	double const pt_btagging=15;

	// Get the jet constituents
	const std::vector<fastjet::PseudoJet>& jet_constituents = jet.constituents();

	// Loop over constituents and look for b quarks
	// also b quarks must be above some minimum pt
	int bquarks = 0;
	for(size_t i=0; i<jet_constituents.size(); i++)
	{
		// Flavour of jet constituent
		const int userid= jet_constituents.at(i).user_index();
		const double pt_bcandidate = jet_constituents.at(i).pt();

		if(abs(userid) ==5 )
			if( pt_bcandidate > pt_btagging)
		  		bquarks++;
	}

 	return bquarks; // no b-jets found
}
