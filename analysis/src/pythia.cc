#include "Pythia8/Pythia.h"
#include "pythia.h"

/*
This routine initializes Pythia8
with all the settings for the shower and underlying event
 */
void InitPythia(Pythia8::Pythia & pythiaRun, string const& eventfile, int const& nevt_max, double& weight_norm )
{
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

  // Verify pythia event sample
  if (pythiaRun.info.nProcessesLHEF() != 1)
  {
    std::cerr << "Error: number of subprocesses in LHE sample is not equal to 1 " <<std::endl;
    std::cerr << "This code does not support this at the moment"<<std::endl;
    exit(-1);
  }

  // Event weight information from Pythia
  const double pythia_wgt = pythiaRun.info.sigmaLHEF(0); // Total sample weight
  weight_norm = pythia_wgt/((double) nevt_max);      // Unit weight

  return;  
}


/*
  Get the information on all final state particles
 */
void get_final_state_particles(Pythia8::Pythia & pythiaRun, finalState& particles, double& unit_weight)
{
  unit_weight = 1;

  // Fetch next record
  if (!pythiaRun.next()) 
    if (pythiaRun.info.atEndOfFile()) // Stop showering when the end of the LHE file is reached
    {
      cerr << "Info: end of input lhe file reached" << endl;
      exit(-1);
    }

  for (int i = 0; i < pythiaRun.event.size(); i++){
    
    // Initialization
    double px = 0;
    double py = 0;
    double pz =0;
    double E = 0;
    
    // Get PDG ID
    int particle_id = pythiaRun.event[i].id();
    
    // Get PDG chargeID (three times the real charge)
   // int particle_charge = pythiaRun.event[i].chargeType();
    
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


