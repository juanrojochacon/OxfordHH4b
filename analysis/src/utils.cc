#include "utils.h"

// General settings for analysis and jet reconstruction
#include "settings.h"

#include "Pythia8/Pythia.h"


/*
This routine initializes Pythia8
with all the settings for the shower and underlying event
 */
void InitPythia(Pythia8::Pythia & pythiaRun, string eventfile){

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
  pythiaRun.readString("PartonLevel:MPI = off"); // Off multiple interactions (UE) 
 
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


/*
This routine initializes Pythia8
with all the settings for the shower and underlying event
 */
void InitPythia_PartonLevel(Pythia8::Pythia & pythiaRun, string eventfile){

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
 
  // No shower
  pythiaRun.readString("SpaceShower:QEDshowerByQ  = off"); // QED shower off
  pythiaRun.readString("SpaceShower:QEDshowerByL  = off"); // QED shower off
  pythiaRun.readString("TimeShower:QEDshowerByQ = off");  // QED off on ISR / quarks irradiate photons
  pythiaRun.readString("TimeShower:QEDshowerByL = off");  // QED off on ISR / leptons irradiate photons  
  pythiaRun.readString("TimeShower:QEDshowerByGamma = off");  // Allow photons to branch into lepton or quark pairs 

  
  // Initial and final state radiation deactivated
  pythiaRun.readString("PartonLevel:ISR = off");  // Shower on
  pythiaRun.readString("PartonLevel:FSR = off");  // Shower on
  
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


/*
  Get the information on all final state particles
 */
void get_final_state_particles(Pythia8::Pythia & pythiaRun, finalState& particles){

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
  
  // Verify final state four momenta
  Analysis::VerifyFourMomentum(particles);
}


// ----------------------------------------------------------------------------------
// Recluster with kt algorithm to obtain splitting scales
std::vector< double > SplittingScales( std::vector<fastjet::PseudoJet> const& jetVec )
{
   
   //vectors that contain the respective splitting scales for all jets
   std::vector<double> split12_vec;
   
   for( int i = 0; i < (int) jetVec.size(); i++){
   
      // For now: Calculate substructure information only for the two leading jets
      if( i > 1 ) continue;
   
      double split12 = -1.;

      if (!jetVec.at(i).has_constituents()){
         std::cout << "ERROR! Splittings (d12, d23, ...) can only be calculated on jets for which the constituents are known."<< std::endl;
         split12_vec.push_back(-1);
         continue;
      }
   
      vector<fastjet::PseudoJet> constits = jetVec.at(i).constituents();
      
//       std::cout << "Jet " <<  i << " has " << constits.size() << " constituents." << std::endl;
      
      fastjet::JetDefinition ekt_jd = fastjet::JetDefinition( fastjet::kt_algorithm, 1.5, fastjet::E_scheme, fastjet::Best);
      const fastjet::ClusterSequence kt_seq_excl = fastjet::ClusterSequence( constits, ekt_jd);
      fastjet::PseudoJet kt_jet = sorted_by_pt( kt_seq_excl.inclusive_jets())[0];
      
      split12 = 1.5*sqrt( kt_seq_excl.exclusive_subdmerge( kt_jet, 1));
      
      split12_vec.push_back(split12);
   }
   
   return split12_vec;
}

// ----------------------------------------------------------------------------------
// Recluster with kt algorithm to obtain nsubjettiness
std::vector< double > NSubjettiness( std::vector<fastjet::PseudoJet> const& jetVec, double const& jet_Rmax, double const& jet_Rmin, double const& jet_Rho )
{

   // Reclustering with VR kt algorithm to obtain nsubjettiness
   
   //vector that contain the respective nsubjettiness variables for all jets
   std::vector<double> tau1_vec;
   std::vector<double> tau2_vec;
   std::vector<double> tau3_vec;
   
   std::vector<double> tau21_vec;
   
   double alpha=1;
   
   for( int i = 0; i < (int) jetVec.size(); i++){
   
      // For now: Calculate substructure information only for the two leading jets
      if( i > 1 ) continue;
   
      double tau1 = -1.;
      double tau2 = -1.;
      double tau3 = -1.;
      
      // Calculate effective jet radius
      double jet_Pt = jetVec[i].pt();
      
      double jet_rad;
      if( jet_Pt > jet_Rmax ) jet_rad = jet_Rmax;
      else if( jet_Pt < jet_Rmin ) jet_rad = jet_Rmin;
      else jet_rad = jet_Rho / jet_Pt;

      if (!jetVec.at(i).has_constituents()){
         std::cout << "ERROR! NSubjettiness (tau1, tau2, ...) can only be calculated on jets for which the constituents are known."<< std::endl;
         
         tau1_vec.push_back(-1);
   tau2_vec.push_back(-1);
   tau3_vec.push_back(-1);
  
   double tau21 = tau2/tau1;
   tau21_vec.push_back(-1);
         continue;
      }
      
      //Code snippet taken from https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetSubstructureVariables#N_subjettiness
      vector<fastjet::PseudoJet> constits = jetVec.at(i).constituents();
      if(constits.size()>0){
      
   double R_kt = 1.5;
   fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::kt_algorithm,R_kt,fastjet::E_scheme,fastjet::Best);
   fastjet::ClusterSequence kt_clust_seq(constits, jet_def);
   vector<fastjet::PseudoJet> kt1axes = kt_clust_seq.exclusive_jets(1);
   double tauNum = 0.0;
   double tauDen = 0.0;
   for (int i = 0; i < (int)constits.size(); i++) {
      // find minimum distance
      double minR = 10000.0; // large number
      for (int j = 0; j < (int)kt1axes.size(); j++) {
         double tempR = sqrt(constits[i].squared_distance(kt1axes[j])); // delta R distance
         if (tempR < minR) minR = tempR;
      }
      tauNum += constits[i].perp() * pow(minR,alpha);
      tauDen += constits[i].perp() * pow(jet_rad,alpha);
   }
   tau1 = tauNum/tauDen;

   if(constits.size()>1){
      vector<fastjet::PseudoJet> kt2axes = kt_clust_seq.exclusive_jets(2);
      tauNum = 0.0;
      tauDen = 0.0;
      for (int i = 0; i < (int)constits.size(); i++) {
      // find minimum distance
      double minR = 10000.0; // large number
      for (int j = 0; j < (int)kt2axes.size(); j++) {
         double tempR = sqrt(constits[i].squared_distance(kt2axes[j]));
         if (tempR < minR) minR = tempR;
      }
      tauNum += constits[i].perp() * pow(minR,alpha);
      tauDen += constits[i].perp() * pow(jet_rad,alpha);
      }
      tau2 = tauNum/tauDen;
      
      if(constits.size() > 2){
      
         vector<fastjet::PseudoJet> kt3axes = kt_clust_seq.exclusive_jets(3);
         tauNum = 0.0;
         tauDen = 0.0;
         for (int i = 0; i < (int)constits.size(); i++) {
         // find minimum distance
         double minR = 10000.0; // large number
         for (int j = 0; j < (int)kt3axes.size(); j++) {
      double tempR = sqrt(constits[i].squared_distance(kt3axes[j]));
      if (tempR < minR) minR = tempR;
         }
         tauNum += constits[i].perp() * pow(minR,alpha);
         tauDen += constits[i].perp() * pow(jet_rad,alpha);
         }
         tau3 = tauNum/tauDen;
      }

   }
      }
      
      tau1_vec.push_back(tau1);
      tau2_vec.push_back(tau2);
      tau3_vec.push_back(tau3);
      
      double tau21 = tau2/tau1;
      tau21_vec.push_back(tau21);
   }
   
   return tau21_vec;
}

