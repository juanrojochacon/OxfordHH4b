/*
This routine simulates the b tagging
In its present form, it only requires that at least
one of the jet constitutents is a b jet
b quarks must be above some minimum pt
 */

bool btagging(int ijet, vector<fastjet::PseudoJet> jets){

  bool btagging_result = false; // Initialization

  // Get the jet constituents
  vector<fastjet::PseudoJet> jet_constituents = jets.at(ijet).constituents();

  // number of b quarks
  int nb=0;

  // Loop over constituents and look for b quarks
  // also b quarks must be above some minimum pt
  for(unsigned i=0;i<jet_constituents.size(); i++){
    int userid= jet_constituents.at(i).user_index();
    double pt_bcandidate = jet_constituents.at(i).pt();
    if(abs(userid) ==5 ){
      // Cut in pt
      if( pt_bcandidate > pt_btagging){
	nb++;
      }
    }
  }

  if(nb > 0 ) btagging_result=true;

  return btagging_result;

}

//-----------------------------------------------------------------

