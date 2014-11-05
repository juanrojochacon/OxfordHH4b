// basic.cc

#include "variableR.h"
#include "utils.h"
#include "settings.h"

// Some useful fastjet includes
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

#include "YODA/Histo1D.h"

using namespace fastjet::contrib;

//----------------------------------------------------------------------
/// a function that pretty prints a list of jets
void print_jets (const fastjet::ClusterSequence & clust_seq,
                 const vector<fastjet::PseudoJet> & jets) {
   
   // sort jets into increasing pt
   vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);
   
   // label the columns
   printf("%5s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity",
          "phi", "pt","m","e", "n constituents");
   
   // print out the details for each jet
   for (unsigned int i = 0; i < sorted_jets.size(); i++) {
      int n_constituents = clust_seq.constituents(sorted_jets[i]).size();
      printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f %8u\n",
             i, sorted_jets[i].rap(), sorted_jets[i].phi(),
             sorted_jets[i].perp(),sorted_jets[i].m(),sorted_jets[i].e(), n_constituents);
   }
}
	
VariableRAnalysis::VariableRAnalysis(string const& sample):
Analysis("variabler", sample)
{
	// Example histo parameters
	const int nbins_example = 20;
	const double histo_min = 0;
	const double histo_max = 500;

	// Example histogram
	BookHistogram(new YODA::Histo1D(nbins_example, histo_min, histo_max), "example_histo");

	// Output nTuple kinematics description (# signal source) is required
	const std::string tupleSpec = "# signal source";
 	outputNTuple<<tupleSpec<<std::endl;	// DO NOT CHANGE
}

void VariableRAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& event)
{
	// Weightnorm provides the sample's unit weight
	double event_weight = weightnorm;

	// Preclustering
	const double rho(2000.), min_r(0.4), max_r(2.0),ptmin(5.0);

	VariableRPlugin lvjet_pluginAKT_precluster(rho, min_r, max_r, VariableRPlugin::AKTLIKE,true);
	fastjet::JetDefinition jet_defAKT_precluster(&lvjet_pluginAKT_precluster);
	fastjet::ClusterSequence clust_seqAKT_precluster(event, jet_defAKT_precluster);

	// tell the user what was done
	cout << "Ran " << jet_defAKT_precluster.description() << endl;

	// extract the inclusive jets with pt > 5 GeV
	vector<fastjet::PseudoJet> inclusive_jetsAKT_precluster = clust_seqAKT_precluster.inclusive_jets(ptmin);

	// print them out
	cout << "Printing inclusive jets with pt > "<< ptmin <<" GeV\n";
	cout << "---------------------------------------\n";
	print_jets(clust_seqAKT_precluster, inclusive_jetsAKT_precluster);
	cout << endl;

	// Write kinematics to nTuple
	outputNTuple <<signal <<"\t"<< GetSample() <<"\t" // DO NOT CHANGE
	<<std::endl; 

	// Pass event
	Pass(event_weight);
}