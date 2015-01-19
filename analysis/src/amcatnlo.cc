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
	Cut("Basic: Two dijets",  0);
	Cut("Leading jet pT", 	  0);
	Cut("Subleading Jet pT",  0);
	Cut("Jet eta acceptance", 0);

	BookHistogram(new YODA::Histo1D(1, 0, 1), "xSec_postCut");

	const std::string tupleSpec = "# signal source";
	outputNTuple<<tupleSpec<<std::endl;
}

void AMCAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
	// Set initial weight
	double event_weight = weightnorm;

	static double const jetR=0.5;
	fastjet::JetDefinition akt(fastjet::antikt_algorithm, jetR);
	fastjet::ClusterSequence cs_akt(fs, akt);

	std::vector<fastjet::PseudoJet> jets_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
	if (!VerifyFourMomentum(jets_akt)) return Cut("4Momentum",event_weight); // Verify clustering

	// We require at least 4 jets in the event, else discard event
	int const njet=4;
	if ( (int)jets_akt.size() < njet) 
		return Cut("Basic: Two dijets",event_weight);

	// Require at least one jet with pT > 100 GeV
	if ( jets_akt[0].pt() < 100 )
		return Cut("Leading jet pT", event_weight);

	// pT / eta cuts
	const double pt_amc_cut = 80;
	const double eta_amc_cut = 2.5;
	for(int ijet=0; ijet<4;ijet++)
	{
		if ( jets_akt[ijet].pt() < pt_amc_cut )
			return Cut("Subleading Jet pT", event_weight);

		if ( fabs(jets_akt[ijet].eta()) > eta_amc_cut )
			return Cut("Jet eta acceptance", event_weight);
	}

	FillHistogram("xSec_postCut", event_weight, 0.5 );

	// ************************************* MVA Output **********************************************************

	outputNTuple <<signal <<"\t"<<GetSample()<<std::endl; 

	// Pass event
	Pass(event_weight);
}


