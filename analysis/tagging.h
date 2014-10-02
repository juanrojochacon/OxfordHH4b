//******************************************************************
//******************************************************************

//-----------------------------------------------------------------

//////////////////////////////
//
// btagging

// To be able to tag a light quark or a b quark
// its pt must the larger than some threshold
//
// We always return four b quarks
// We give a weigh for each event

bool btagging(int ijet, vector<fastjet::PseudoJet> jets){

  bool btagging_result = false; // Initialization
  
  // Get the jet constituents
  vector<fastjet::PseudoJet> jet_constituents = jets.at(ijet).constituents();

  // number of b quarks
  int nb=0;

  // Loop over constituents and look for b quarks
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

  for(int ib=0; ib< nb; ib++){
     
    // Now look over number of b quarks
    double random = double(rand())/RAND_MAX;
    if(random < 0.0 || random > 1.0){
      std::cout<<"Invalid value of random number"<<std::endl;
      exit(-10);
    }
    
    if(random < btag_prob) btagging_result=true;
    if(btagging_result) break;

  }

  // Here mistag probability, in the case there is no b quarks
  if(nb==0 ){

    // Check that at least one light parton above some pt threshold
    bool  pt_bcandidate_check=false;
    for(unsigned i=0;i<jet_constituents.size(); i++){
      double pt_bcandidate = jet_constituents.at(i).pt();
      if( pt_bcandidate > pt_btagging)  pt_bcandidate_check=true;
    }
    if(pt_bcandidate_check) btagging_result=false;
    return btagging_result;

    double random = double(rand())/RAND_MAX;
    if(random < 0.0 || random > 1.0){
      std::cout<<"Invalid value of random number"<<std::endl;
      exit(-10);
    }

    // mistag probabily of quark and gluon jets
    if(random < btag_mistag) {
      btagging_result=true;
      // std::cout<<"light jet mistag"<<std::endl;
      //std::cout<< jets.at(ijet).pt()<<std::endl;
    }
  }
  
  // Return the result of the b tagging
  // std::cout<<btagging_result<<std::endl;
  return btagging_result;

}



//-----------------------------------------------------------------
//-----------------------------------------------------------------
