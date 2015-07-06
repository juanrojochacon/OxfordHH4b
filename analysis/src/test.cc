// oxford_res_fr.cc

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

#include "testCluster.h"
#include "utils.h"
#include "settings.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

#include <algorithm>

using namespace fastjet::contrib;

// Resolved jet radius
const double ResJetR=0.4;


TestAnalysis::TestAnalysis(std::string const& sampleName):
Analysis("test", sampleName)
{
  
  // ********************* Ntuple definition **********************
  const std::string tupleSpec = "# signal source weight pt_H0 pt_H1 pt_HH m_H0 m_H1 m_HH dR_HH dPhi_HH dEta_HH";
  outputNTuple<<tupleSpec<<std::endl;
}

void OxfordCombinedRWAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& fs)
{
  Analysis::Analyse(signal, weightnorm, fs);

  // Set initial weight
  const double event_weight = weightnorm;

    // ********************************* Simple event categorisation  ***************************************

  fastjet::JetDefinition akt_res(fastjet::antikt_algorithm, ResJetR);
  fastjet::ClusterSequence cs_akt_res(fs, akt_res);
  std::vector<fastjet::PseudoJet> smallRJets = sorted_by_pt( cs_akt_res.inclusive_jets()  ); // Get all the jets (no pt cut here)
  
  // Basic kinematic cuts
  std::vector<fastjet::PseudoJet> smallRJetsSel;
  std::vector<fastjet::PseudoJet> smallRJetsSel_unsrt;
  for( size_t i = 0; i < smallRJets.size(); i++)
  {
    cout << GetSample() << "  "<< smallRJets.at(i).pt() <<"  "<< fabs( smallRJets.at(i).eta() ) <<endl;
    if( smallRJets.at(i).pt() < 40. ) continue;
    if( fabs( smallRJets.at(i).eta() ) > 2.5 ) continue;
    
    smallRJetsSel_unsrt.push_back( smallRJets.at(i) );
  }
 }