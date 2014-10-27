// newanalyis.h jr/nh 21/10/14
#pragma once

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>

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
		Analysis(string const& name, string const& sample);
		virtual ~Analysis();

		string const& GetName() const {return analysisName;};
		string const& GetRoot() const {return analysisRoot;};
		string const& GetSample() const {return sampleName;};

		int const& GetNPassed() const {return nPassed;};
		double const& GetWeight() const {return passedWeight;};

		virtual void Analyse(bool const& signal, double const& weightnorm, finalState const&) = 0;
		void Export();

	protected:
		void BookHistogram(YODA::Histo1D*, string const& name);
		void FillHistogram(string const& rname, double const& weight, double const& coord );

		void Cut(string const& cutname, double const& weight);
		void Pass(double const& weight);

		std::ofstream outputNTuple;	//!< Ntuple record for all events

	private:
		// Standard analysis info
		const string analysisName;		//!< Name of analysis
		const string analysisRoot;		//!< Path to root of analysis results folder
		const string sampleName;		//!< Name of the current sample

		int nPassed; 			//!< Number of events passing analysis cuts
		double passedWeight; 	//!< Total weight of passed events

		std::map<int,YODA::Histo1D*> bookedHistograms;
		std::vector<std::pair<const std::string, double> > cutWeight;
};

