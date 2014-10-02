//******************************************************************
//******************************************************************

bool bbbb_analysis(vector<fastjet::PseudoJet> jets_akt){
  
  bool hhtag=true;


  int const njet=4;
  // Exit if less than four jets
  if(jets_akt.size() < njet) return false;

  // Identify the four b jets
  // if failed, returns false
  
  vector<fastjet::PseudoJet> bjets;
  for(unsigned ijet=0; ijet<jets_akt.size();ijet++){
  
    if(jets_akt.at(ijet).pt() < pt_btagging) continue;
    bool btag = btagging(ijet, jets_akt);

    if(btag) bjets.push_back(jets_akt.at(ijet));

  }
  
  // Exit if less that 4 b jets found
  bjets=sorted_by_pt(bjets);
  if(bjets.size() < njet) return false;
 
  // Check that the four leading b jets are in acceptance
  for(unsigned ijet=0; ijet<njet;ijet++){

    if(bjets.at(ijet).pt() < pt_bjet) return false;
    if(fabs( bjets.at(ijet).eta() ) > eta_bjet) return false;
    
  }

  // Now pair them to produce the Higgs mass
  // Get the pairing that minimizes |m_dj1 - m_dj2|
  double dijet_mass[njet][njet];
  for(unsigned ijet=0;ijet<njet;ijet++){
    for(unsigned jjet=0;jjet<njet;jjet++){
      // Compute jet masses
      PseudoJet sum;
      sum = bjets[ijet] + bjets[jjet];
      dijet_mass[ijet][jjet] = sum.m();
    }
  }
	      
  double mdj_diff_min = 1e20; // Some large number to begin
  unsigned jet1_id1=10,jet1_id2=10,jet2_id1=10,jet2_id2=10;
	      
  for(unsigned ijet=0;ijet<njet;ijet++){
    for(unsigned jjet=ijet+1;jjet<njet;jjet++){
      double mdj1 = dijet_mass[ijet][jjet];
      for(unsigned ijet2=0;ijet2<njet;ijet2++){
	for(unsigned jjet2=ijet2+1;jjet2<njet;jjet2++){
	  double mdj2 = dijet_mass[ijet2][jjet2];
	  double min_dj = fabs(mdj1 - mdj2);
	  if(min_dj <  mdj_diff_min && ijet != ijet2  && jjet != jjet2 && jjet !=ijet2 && ijet != jjet2 ){
	    mdj_diff_min = min_dj;
	    jet1_id1 = ijet;
	    jet1_id2 = jjet;
	    jet2_id1 = ijet2;
	    jet2_id2 = jjet2;
	  }
	}
      }
    }
  }
  
  // Construct the Higgs candidates
  PseudoJet higgs1 = bjets.at( jet1_id1) + bjets.at( jet1_id2); 
  PseudoJet higgs2 = bjets.at( jet2_id1) + bjets.at( jet2_id2);
	    
  // Higgs mass window condition
  double mass_diff1 = fabs(higgs1.m()-m_higgs)/m_higgs;
  double mass_diff2 = fabs(higgs2.m()-m_higgs)/m_higgs;

  if( mass_diff1 > mass_resolution || mass_diff2 > mass_resolution ) {
    return false;
  }

  //std::cout<<higgs1.m()<<" "<<higgs2.m()<<std::endl;

  // Here fill the histogram for the pt of the hh system
  string histofill="pthh";
  double xsec_fill=1.0;
  PseudoJet dihiggs = higgs1+higgs2;
  double pthh = dihiggs.pt();
  histo_fill(histofill, xsec_fill, pthh);

  // Fill histogram for pt of individual higgses
  histofill="pth";
  double pth = higgs1.pt();
  histo_fill(histofill, xsec_fill, pth);
  pth = higgs2.pt();
  histo_fill(histofill, xsec_fill, pth);

  // Now impose kinematical cut in the pt of the hh system
  if(pthh < pthh_cut) return false;

  return hhtag;

}



//******************************************************************
//******************************************************************

bool bbbb_analysis_boosted(vector<fastjet::PseudoJet> jets_akt){
  
  bool hhtag=true;


  int const njet=2;
  // Exit if less than two jets
  // We are in the fully boosted regime
  if(jets_akt.size() < njet) return false;

  // Now look for substructure in the two hardest jets
  vector<PseudoJet> higgs;

  for (unsigned i = 0; i < njet; i++) {
    
    // first recluster with some large CA (needed for mass-drop)
    ClusterSequence cs_sub(jets_akt[i].constituents(), CA10);
    
    // next get hardest jet
    PseudoJet ca_jet = sorted_by_pt(cs_sub.inclusive_jets())[0];
    
    // now run mass drop tagger
    MassDropTagger md_tagger(mu, ycut);
    PseudoJet tagged_jet = md_tagger(ca_jet);
    
    // Run also b-tagging
    bool btagging_result = btagging(i,jets_akt);
    
    if (tagged_jet == 0 ) {
	    
      // No mass-drop tag
      
    } else {
	    
      // Assign as a jet if on top of MDT also b tagging is positive
      // Here we would need to implement double b-tagging to 
      // improve background rejection
      if(btagging_result) higgs.push_back(tagged_jet);
	      
    }

  }

  // If only and two only jets have been both mass-drop tagged
  // and b-tagged, continute
  if(higgs.size()!=2) return false;

   // Higgs mass window condition
  double mass_diff1 = fabs(higgs.at(0).m()-m_higgs)/m_higgs;
  double mass_diff2 = fabs(higgs.at(1).m()-m_higgs)/m_higgs;

  if( mass_diff1 > mass_resolution || mass_diff2 > mass_resolution ) {
    return false;
  }

  //std::cout<<higgs1.m()<<" "<<higgs2.m()<<std::endl;

  // Here fill the histogram for the pt of the hh system
  string histofill="pthh";
  double xsec_fill=1.0;
  PseudoJet dihiggs = higgs.at(0)+higgs.at(1);
  double pthh = dihiggs.pt();
  histo_fill(histofill, xsec_fill, pthh);

  // Fill histogram for pt of individual higgses
  histofill="pth";
  double pth = higgs.at(0).pt();
  histo_fill(histofill, xsec_fill, pth);
  pth = higgs.at(1).pt();
  histo_fill(histofill, xsec_fill, pth);

  // Now impose kinematical cut in the pt of the hh system
  if(pthh < pthh_cut) return false;

  return hhtag;

}
