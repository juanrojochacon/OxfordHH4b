// amcatnlo.cc

#include "amcatnlo.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"

#include "YODA/Histo1D.h"

AMCAnalysis::AMCAnalysis(std::string const& sampleName):
Analysis("aMC@NLO", sampleName)
{
	// CutFlow
	Cut("Basic: Four jets",  0);
	Cut("Basic: Leading jet pT", 	  0);
	Cut("Unrec Sample", 0);
	Cut("NbJets 80",  0);
	Cut("NbJets 100", 0);

	BookHistogram(new YODA::Histo1D(1, 0, 1), "xSec_postCut");

	const std::string tupleSpec = "# signal source weight";
	outputNTuple<<tupleSpec<<std::endl;
}

void AMCAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
	Analysis::Analyse(signal, weightnorm, fs);

	// Set initial weight
	double event_weight = weightnorm;

	const double jetR = 0.5;
	const double eta_amc_cut = 2.5;

	fastjet::JetDefinition akt(fastjet::antikt_algorithm, jetR);
	fastjet::ClusterSequence cs_akt(fs, akt);

	std::vector<fastjet::PseudoJet> jets_akt_80  = sorted_by_pt( cs_akt.inclusive_jets( 80.0  )  );
	std::vector<fastjet::PseudoJet> jets_akt_100 = sorted_by_pt( cs_akt.inclusive_jets( 100.0 )  );

	// Need at least 4 jets with pT > 80
	if ( (int)jets_akt_80.size() < 4) 
		return Cut("Basic: Four jets",event_weight);

	// Need at least 1 jet with pT > 100
	if ( (int)jets_akt_100.size() < 1) 
		return Cut("Basic: Leading jet pT",event_weight);

	// b-Tagging
	int NbJet_80 = 0;
	int NbJet_100 = 0;

	for (size_t i=0; i<jets_akt_80.size(); i++)
		if ( fabs(jets_akt_80[i].eta()) < eta_amc_cut )
			if (BTagging(jets_akt_80[i]))
				NbJet_80++;

	for (size_t i=0; i<jets_akt_100.size(); i++)
		if ( fabs(jets_akt_100[i].eta()) < eta_amc_cut )
			if (BTagging(jets_akt_100[i]))
				NbJet_100++;

	int Nb_Req_80 = 0;
	int Nb_Req_100 = 0;

	if (GetSample().compare("SHERPA_QCD4b") == 0)
	{
		Nb_Req_80 = 4;
		Nb_Req_100 = 1;
	} else if (GetSample().compare("SHERPA_QCD2b2j") == 0)
	{
		Nb_Req_80 = 2;
		Nb_Req_100 = 1;
	}
	else if (GetSample().compare("SHERPA_QCD4j") == 0)
	{
		Nb_Req_80 = 0;
		Nb_Req_100 = 0;
	} else
	{
		return Cut("Unrec Sample", event_weight);
	}

	// Require either 0, 2 or 4 b-Jets at 80GeV
	if ( NbJet_80 < Nb_Req_80 )
		return Cut("NbJets 80", event_weight);

	// Require wither 1 or 0 b-Jets at 100GeV
	if ( NbJet_100 < Nb_Req_100 )
		return Cut("NbJets 100", event_weight);


	FillHistogram("xSec_postCut", event_weight, 0.5 );

	// ************************************* MVA Output **********************************************************

	outputNTuple <<signal <<"\t"<< GetSample() <<"\t"<<event_weight<<"\t"<<std::endl; 

	// Pass event
	Pass(event_weight);
}

bool AMCAnalysis::BTagging( fastjet::PseudoJet const& jet ) const
{
	// Get the jet constituents
	const std::vector<fastjet::PseudoJet>& jet_constituents = jet.constituents();

	// Loop over constituents and look for b quarks
	// also b quarks must be above some minimum pt
	for(size_t i=0; i<jet_constituents.size(); i++)
	{
		// Flavour of jet constituent
		const int userid= jet_constituents.at(i).user_index();
		if(abs(userid) ==5 )
	  		return true;
	}

 	return false; 
}


