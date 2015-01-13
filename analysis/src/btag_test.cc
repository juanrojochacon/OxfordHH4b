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

bTagTestAnalysis::bTagTestAnalysis(std::string const& sampleName):
Analysis("bTagTest", sampleName)
{
	const std::string tupleSpec = "# signal source";
	outputNTuple<<tupleSpec<<std::endl;

	// Order cutflow
	Cut("Basic: Two dijets", 0);
	Cut("Basic: bTagging", 0);
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
	for(int ijet=0; ijet<njet;ijet++)
		if( BTagging(jets_fr_akt[ijet]) )   // Check if at least one of its constituents is a b quark
		{
			bjets.push_back(jets_fr_akt.at(ijet));
			event_weight *= test_btag_prob; // Account for b tagging efficiency
		}
		else // Else, account for the fake b-tag probabililty
		{
			bjets.push_back(jets_fr_akt.at(ijet));
			event_weight *= test_btag_mistag;
		}

	// cut from btagging
	Cut("Basic: bTagging", initial_weight - event_weight);

	
	// ************************************* MVA Output **********************************************************
	outputNTuple <<signal <<"\t"<<GetSample()<<std::endl; 

	// Pass event
	Pass(event_weight);

}

// ----------------------------------------------------------------------------------

bool bTagTestAnalysis::BTagging( fastjet::PseudoJet const& jet ) const
{
	// Cuts for the b-jet candidates for b-tagging
	double const pt_btagging=0;

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


