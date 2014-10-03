/*
This routine initializes Pythia8
with all the settings for the shower and underlying event
 */
void InitPythia(Pythia & pythiaRun, string eventfile){

  // Initialize random seed
  srand (time(NULL));
  std::cout<<"time = "<<time(NULL)<<std::endl;
  double random = double(rand())/RAND_MAX;
  std::cout<<"\n\n Random number I = "<<random<<"\n\n"<<std::endl;
  random = double(rand())/RAND_MAX;
  std::cout<<"\n\n Random number II = "<<random<<"\n\n"<<std::endl;
  
  // Random seed
  pythiaRun.readString("Random:setSeed = on");
  double random_seed_pythia = 100000 * double(rand())/RAND_MAX;
  ostringstream o;
  o<<"Random:seed = "<<int(random_seed_pythia);
  cout<<o.str()<<endl;
  pythiaRun.readString(o.str());

  // Initialize Les Houches Event File run. List initialization information.
  pythiaRun.readString("Beams:frameType = 4"); 
  
  // Shower settings
 
  // The shower is QCD only, no QED or weak effects included
  pythiaRun.readString("SpaceShower:QEDshowerByQ  = off"); // QED shower off
  pythiaRun.readString("SpaceShower:QEDshowerByL  = off"); // QED shower off
  pythiaRun.readString("TimeShower:QEDshowerByQ = off");  // QED off on ISR / quarks irradiate photons
  pythiaRun.readString("TimeShower:QEDshowerByL = off");  // QED off on ISR / leptons irradiate photons  
  pythiaRun.readString("TimeShower:QEDshowerByGamma = off");  // Allow photons to branch into lepton or quark pairs 
  // Initial and final state radiation activated
  pythiaRun.readString("PartonLevel:ISR = on");  // Shower on
  pythiaRun.readString("PartonLevel:FSR = on");  // Shower on
  
  // No hadronization
  pythiaRun.readString("HadronLevel:all = off"); // Of hadronization
 
  // For the time being no  UE or PU included
  pythiaRun.readString("PartonLevel:MI = off"); // Off multiple interactions (UE) 
 
  // Higgs decays always into 4b
  // Need to correct by hand the xsecs for the BR(HH->4b) branching fraction
  pythiaRun.readString("25:onMode = off");
  pythiaRun.readString("25:onIfAll = 5 -5");

  // b quarks and do not decay
  // They are treated as stable particles in the detector
  pythiaRun.readString("5:mayDecay = no");
  pythiaRun.readString("-5:mayDecay = no");

  // Read the Les Houches Event File
  string ofile;
  ofile="Beams:LHEF = "+eventfile;
  pythiaRun.readString(ofile.c_str());

   // Main initialization
  pythiaRun.init();

  std::cout<<"\n Pythia8 Initialized \n "<<std::endl;
  
}

//-------------------------------------------------------------------

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

  // b-tagging is succesful if at least one b quark with
  // pt > pt_btagging has been found
  if(nb > 0 ) btagging_result=true;

  return btagging_result;

}

//----------------------------------------------------------------------------

/*
  Get the information on all final state particles
 */

void get_final_state_particles(Pythia & pythiaRun, vector<fastjet::PseudoJet> & particles){

  for (int i = 0; i < pythiaRun.event.size(); i++){
    
    // Initialization
    double px = 0;
    double py = 0;
    double pz =0;
    double E = 0;
    
    // Get PDG ID
    int particle_id = pythiaRun.event[i].id();
    
    // Get particle status: in pythia8, status > 0 means final state particles
    int particle_status = pythiaRun.event[i].status();
    
    // Consider only final state particles
    if( particle_status <= 0 ) continue;
    
    // Get the particle kinematics
    px= pythiaRun.event[i].px();
    py= pythiaRun.event[i].py();
    pz= pythiaRun.event[i].pz();
    E= pythiaRun.event[i].e();
    
    // quarks and gluons
    // including b-quarks
    if(abs(particle_id)<6 || particle_id==21 ){
      // Set kinematics
      particles.push_back( fastjet::PseudoJet(px,py,pz,E) );
      // Set PDG ID
      particles.at(particles.size()-1).set_user_index(particle_id);
    }
    // beam remnants
    else if(particle_id > 2000 ){
      particles.push_back( fastjet::PseudoJet(px,py,pz,E) );
      particles.at(particles.size()-1).set_user_index(particle_id);
    }
    else{
      std::cout<<"Invalid particle ID = "<<particle_id<<std::endl;
      exit(-10);
    }
    
  } // End loop over particles in event
  
  // Check event energy conservation here
  double px_tot=0;
  double py_tot=0;
  double pz_tot=0;
  double E_tot=0;
  
  // Loop over all particles
  for(unsigned ij=0;ij<particles.size();ij++){
    px_tot+= particles.at(ij).px();
    py_tot+= particles.at(ij).py();
    pz_tot+= particles.at(ij).pz();
    E_tot+= particles.at(ij).E();
  }
  // Check energy-momentum conservation
  if( fabs(px_tot) > tol_emom || fabs(py_tot)  > tol_emom 
      || fabs(pz_tot)  > tol_emom || fabs(E_tot-Eref)  > tol_emom ){
    std::cout<<"\n ********************************************************************** \n"<<std::endl;
    std::cout<<"No conservation of energy in Pythia after shower "<<std::endl;
    std::cout<<"px_tot = "<<px_tot<<std::endl;
    std::cout<<"py_tot = "<<py_tot<<std::endl;
    std::cout<<"pz_tot = "<<pz_tot<<std::endl;
    std::cout<<"E_tot, Eref = "<<E_tot<<" "<<Eref<<std::endl;
    exit(-10);
    std::cout<<"\n ********************************************************************** \n"<<std::endl;
  }
}


//-------------------------------------------------------------------

/*
This routine read the event kinematics and performs the jet clustering
It also checkes that energy momentum is conserved event by event
This applies for small R jet clustering with the anti-kt algorithm
 */

void jet_clustering_analysis_smallR(Pythia & pythiaRun, vector<fastjet::PseudoJet> & bjets, double & event_weight ){

   // This are all the particles that enter the jet clustering
  vector<fastjet::PseudoJet> particles;
  // Get all final state particles
  get_final_state_particles(pythiaRun, particles);
  
  // Now perform jet clustering with anti-kT
  // jetR is defined in settings.h
  // Note that here we use a small R clustering
  JetDefinition akt(antikt_algorithm, jetR_0p5);
  // Cluster all particles
  // The cluster sequence has to be saved to be used for jet substructure
  ClusterSequence cs_akt(particles, akt);
  // Get all the jets (no pt cut here)
  vector<fastjet::PseudoJet> jets_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
  
  // Check again four-momentum conservation, this time applied to jets
  // formed from the clustering of quarks and gluons (and beam remnants as well)
  double px_tot=0;
  double py_tot=0;
  double pz_tot=0;
  double E_tot=0;
  for(unsigned ij=0;ij<jets_akt.size();ij++){
    px_tot+= jets_akt.at(ij).px();
    py_tot+= jets_akt.at(ij).py();
    pz_tot+= jets_akt.at(ij).pz();
    E_tot+= jets_akt.at(ij).E();
  }
  
  // Check energy-momentum conservation
  if( fabs(px_tot) > tol_emom || fabs(py_tot)  > tol_emom 
      || fabs(pz_tot)  > tol_emom || fabs(E_tot-Eref)  > tol_emom ){
    std::cout<<"\n ********************************************************************** \n"<<std::endl;
    std::cout<<"No conservation of energy in Pythia after shower and jet reconstruction "<<std::endl;
    std::cout<<"px_tot = "<<px_tot<<std::endl;
    std::cout<<"py_tot = "<<py_tot<<std::endl;
    std::cout<<"pz_tot = "<<pz_tot<<std::endl;
    std::cout<<"E_tot, Eref = "<<E_tot<<" "<<Eref<<std::endl;
    exit(-10);
    std::cout<<"\n ********************************************************************** \n"<<std::endl;
  }
  
  // Now Initialize the event weight
  event_weight=1.0;

  // We require at least 4 jets in the event, else discard event
  int const njet=4;
  if(jets_akt.size() < njet) {
    event_weight=0;
    return;
  }

  // By looking at the jet constituents
  // we can simulate the effects of b tagging

  // Loop over the 4 hardest jets in event only
  for(unsigned ijet=0; ijet<njet;ijet++){
  
    // Check if at least one of its constituents is a b quark
    bool btag = btagging(ijet, jets_akt);
    
    // If the jet has a pt > pt_btagging, assign this jet to be
    // a b jet
    if(jets_akt.at(ijet).pt() > pt_btagging){

      // Account for b tagging efficiency
      if(btag) {
	bjets.push_back(jets_akt.at(ijet));
	event_weight *= btag_prob;
      }
      // Else, account for the fake b-tag probabililty
      else{
	bjets.push_back(jets_akt.at(ijet));
	event_weight *= btag_mistag;
      }
      // std::cout<<"ijet, btag = "<<ijet<<"\t"<<btag<<std::endl;
      
    }

  }
    
  // Exit if less that 4 b jets found
  // Recall that b jets require a minimum pt 
  // to simulate realistic b tagging
  bjets=sorted_by_pt(bjets);
  if(bjets.size() < njet) {
    event_weight=0;
    return;
  }


} 

// ----------------------------------------------------------------------------------

/*
This routine read the event kinematics and performs the jet clustering
It also checkes that energy momentum is conserved event by event
This applies for large R R jet clustering with the anti-kt algorithm
Here we tag large-R jets using jet substructure tools 
 */

void jet_clustering_analysis_largeR(Pythia & pythiaRun, vector<fastjet::PseudoJet> & higgs_candidates, double & event_weight ){

  // std::cout<<"in jet_clustering_analysis_largeR"<<std::endl;

   // This are all the particles that enter the jet clustering
  vector<fastjet::PseudoJet> particles;
  // Get all final state particles
  get_final_state_particles(pythiaRun, particles);
  
  // Now perform jet clustering with anti-kT
  // jetR is defined in settings.h
  // Here we use a large-R clustering, R=1.2
  JetDefinition akt(antikt_algorithm, jetR_1p2);
  // Cluster all particles
  // The cluster sequence has to be saved to be used for jet substructure
  ClusterSequence cs_akt(particles, akt);
  // Get all the jets (no pt cut here)
  vector<fastjet::PseudoJet> jets_akt = sorted_by_pt( cs_akt.inclusive_jets()  );
  
  // Check again four-momentum conservation, this time applied to jets
  // formed from the clustering of quarks and gluons (and beam remnants as well)
  double px_tot=0;
  double py_tot=0;
  double pz_tot=0;
  double E_tot=0;
  for(unsigned ij=0;ij<jets_akt.size();ij++){
    px_tot+= jets_akt.at(ij).px();
    py_tot+= jets_akt.at(ij).py();
    pz_tot+= jets_akt.at(ij).pz();
    E_tot+= jets_akt.at(ij).E();
  }
  if( fabs(px_tot) > tol_emom || fabs(py_tot)  > tol_emom 
      || fabs(pz_tot)  > tol_emom || fabs(E_tot-Eref)  > tol_emom ){
    std::cout<<"\n ********************************************************************** \n"<<std::endl;
    std::cout<<"No conservation of energy in Pythia after shower and jet reconstruction "<<std::endl;
    std::cout<<"px_tot = "<<px_tot<<std::endl;
    std::cout<<"py_tot = "<<py_tot<<std::endl;
    std::cout<<"pz_tot = "<<pz_tot<<std::endl;
    std::cout<<"E_tot, Eref = "<<E_tot<<" "<<Eref<<std::endl;
    exit(-10);
    std::cout<<"\n ********************************************************************** \n"<<std::endl;
  }
  
  // Now initialize the event weight
  event_weight=1.0;

  // We require at least 2 large-R jets in the event, else discard event
  int const njet=2;
  if(jets_akt.size() < njet) {
    event_weight=0;
    return;
  }

  // The two hardest jets must have pt > 200 GeV
  // Before doing the cut fill the ptH histogram
  histo_fill("pth", event_weight, jets_akt.at(0).pt());
  histo_fill("pth", event_weight, jets_akt.at(1).pt());

  // Now the cut on the pt of the Higgs candidates
  double const pt_largeRjet = 200.0;
  for(unsigned ijet=0;ijet<njet;ijet++){
    if(jets_akt.at(ijet).pt() < pt_largeRjet) {
      event_weight=0;
      return;
    }
  }

  // Now look for substructure in each of these two dijets using the BDRS
  // mass-drop tagger
  
  for (unsigned i = 0; i < njet; i++) {
    
    // first recluster again with the large-R but with Cambridge-Aachen
    ClusterSequence cs_sub(jets_akt[i].constituents(), CA10);
    // get hardest jet
    PseudoJet ca_jet = sorted_by_pt(cs_sub.inclusive_jets())[0];
    
    // now run mass drop tagger
    // parameters are specified in settings.h
    MassDropTagger md_tagger(mu, ycut);
    PseudoJet tagged_jet = md_tagger(ca_jet);
    
    // If tagging succesful - declare as Higgs candidate
    if (tagged_jet != 0 )  higgs_candidates.push_back(tagged_jet);
    
  }

  // If we don't have a mass-drop tag in each of the two leading large-R jets
  // discard the event
  if(higgs_candidates.size()!=2) {
    event_weight=0.0;
    return;
  }

  // Now we need to implement b-tagging
  // Here different options possible
  // To begin with, most agressive option
  // Require only two b quarks in the fat jet
  // each with pt > 40 GeV
  // No restruction on how close they can be in angle
  double const pt_btagging_largeR=40.0;

  // Get the jet constituents
  int const nb_fatjet=2;
  for (unsigned i = 0; i < njet; i++) {

    vector<fastjet::PseudoJet> jet_constituents = jets_akt.at(i).constituents();

    // number of b quarks
    int nb=0;

    // Loop over constituents and look for b quarks
    // These b quarks must be hard enough
    for(unsigned i=0;i<jet_constituents.size(); i++){
      int userid= jet_constituents.at(i).user_index();
      double pt_bcandidate = jet_constituents.at(i).pt();
      if(abs(userid) ==5 ){
	// Cut in pt
	if( pt_bcandidate > pt_btagging_largeR){
	  nb++;
	}
      }
    }
    // Now assign the event weight
    // This can be refined, for example requiring for fakes to have subjets with
    // pt above some threshold
    if(nb >= nb_fatjet) event_weight *= pow(btag_prob,2.0);
    if(nb ==1 ) event_weight *= (btag_prob*btag_mistag);
    if(nb ==0 ) event_weight *= pow(btag_mistag,2.0);
  }

  //std::cout<<"jet_clustering_analysis_largeR completed"<<std::endl;

  return;

} 

// ----------------------------------------------------------------------------------

/*
This is the analysis used by the UCL group
See for example the slides of their talk at Boost 2014
https://indico.cern.ch/event/302395/session/12/contribution/26/material/slides/1.pdf
 */

bool analysis_4b_ucl(vector<fastjet::PseudoJet> & bjets, double event_weight ){

  // First of all, after basic selection, require that all four b jets are above 40 GeV
  double const pt_bjet_ucl = 40.0;
  // they should also be in central rapodity, |eta| < 2.5
  double const eta_bjet_ucl = 2.5;

  int const njet=4;
  // Restrict to the four leading jets in the event
  for(unsigned ijet=0; ijet<njet;ijet++){
    // std::cout<<"ijet, pt, eta = "<<ijet<<" \t "<<bjets.at(ijet).pt()<<" \t "<<bjets.at(ijet).eta()<<std::endl;
    if(bjets.at(ijet).pt() < pt_bjet_ucl || 
       fabs( bjets.at(ijet).eta() ) > eta_bjet_ucl) return false;
  }
  
  // std::cout<<"4 b jets passing ucl acceptance found"<<std::endl;

  // The next step is to apply the dijet selection
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

  // Histograms before cut in pt Higgs candidates
  histo_fill("pth", event_weight, higgs1.pt());
  histo_fill("pth", event_weight, higgs2.pt());

  // Now, the pt of the dijet Higgs candidates must be above 150 GeV
  double pt_dijet_ucl=150.0;
  if(higgs1.pt() < pt_dijet_ucl ||higgs1.pt() < pt_dijet_ucl ) return false;

  // std::cout<<"Higgs pt cut passed"<<std::endl;

  // These two dijets cannot be too far in DeltaR
  // Check exactly the cut used: deltaR or delta_eta?
  double delta_eta_dijet_ucl=1.5;
  double delta_eta_dijet = fabs(higgs1.eta()- higgs2.eta());
  //cout<<"Delta_eta_dijet = "<<delta_eta_dijet<<std::endl;
  if(delta_eta_dijet > delta_eta_dijet_ucl) return false;

  //std::cout<<"Higgs deltaR cut passed"<<std::endl;
	    
  // Higgs mass window condition
  double mass_diff1 = fabs(higgs1.m()-m_higgs)/m_higgs;
  double mass_diff2 = fabs(higgs2.m()-m_higgs)/m_higgs;
  if( mass_diff1 > mass_resolution || mass_diff2 > mass_resolution ) return false;

  // Histograms for the pt of the HH system
  // no cuts are applied on this variable
  PseudoJet dihiggs= higgs1+higgs2;
  histo_fill("pthh", event_weight, dihiggs.pt());
  
  //  std::cout<<"Event tagged as HH->4b event"<<std::endl;

  return true;

}

//---------------------------------------------------------------------------------


/*
This analysis is inspired by the Durham group, based on jet substructure
See for example the slides of their talk at Boost 2014
https://indico.cern.ch/event/302395/session/12/contribution/25/material/slides/0.pdf
as well as more details in their publication
 arXiv:1404.7139
Here we only use BDRS tagging, no attemp to use Shower Deconstruction
 */

bool analysis_4b_durham(vector<fastjet::PseudoJet> & higgs_candidates, double event_weight ){

  // std::cout<<"in analysis_4b_durham"<<std::endl;

  // Require that the Higgs candidates are in acceptance
  double const eta_bjet_durham = 2.5;
  
  int const njet=2;
  for(unsigned ijet=0; ijet<njet;ijet++){
    if(fabs( higgs_candidates.at(ijet).eta() ) > eta_bjet_durham) return false;
  }

  // Same as in the UCL analysis
  // require that these two leading fat jets are not too separated in rapidity
  double delta_eta_dijet_fatjet=1.5;
  double delta_eta_dijet = fabs(higgs_candidates.at(0).eta()- higgs_candidates.at(1).eta());
  if(delta_eta_dijet > delta_eta_dijet_fatjet) return false;
  
  // Higgs mass window condition
  double mass_diff1 = fabs(higgs_candidates.at(0).m()-m_higgs)/m_higgs;
  double mass_diff2 = fabs(higgs_candidates.at(1).m()-m_higgs)/m_higgs;
  if( mass_diff1 > mass_resolution || mass_diff2 > mass_resolution ) return false;
  
  // Histograms for the pt of the HH system
  // no cuts are applied on this variable
  PseudoJet dihiggs= higgs_candidates.at(0)+higgs_candidates.at(1);
  histo_fill("pthh", event_weight, dihiggs.pt());
  
  //std::cout<<"Event tagged as HH->4b event"<<std::endl;

  return true;

}

//--------------------------------------------------------------------------------------------


