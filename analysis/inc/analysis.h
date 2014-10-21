// newanalyis.h jr/nh 21/10/14
#pragma once

#include <iostream>
#include <fstream>
#include <list>
#include <vector>

using std::string;

// Fwd fastjet
namespace fastjet{
	class PseudoJet;
}

// Fwd YODA
namespace YODA{
	class Histo1D;
}

// Final state particles
typedef std::vector<fastjet::PseudoJet> finalState;

class Analysis
{
	public:
		Analysis(string const& name);
		~Analysis();

		string const& GetName() const {return analysisName;};
		string const& GetRoot() const {return analysisRoot;};

		void InitSample(string const& sampleID);

		virtual void Analyse(string const& sampleID, bool const& signal, finalState const&) = 0;

	protected:
		void BookHistogram(YODA::Histo1D*, string const& name);
		void FillHistogram(string const& rname, string const& sample, double const& weight, double const& coord );

		int nPassed; 			//!< Number of events passing analysis cuts
		double passedWeight; 	//!< Total weight of passed events

		std::ofstream totalNTuple;	//!< Ntuple record for all events
		std::ofstream sampleNTuple;	//!< Ntuple record per LHE sample

	private:
		// Standard analysis info
		const string analysisName;		//!< Name of analysis
		const string analysisRoot;		//!< Path to root of analysis results folder

		std::list<YODA::Histo1D*> protoHistograms; //!< Prototype histograms
		std::list<string> bookedHistograms; 	   //!< Histograms booked into system
};

