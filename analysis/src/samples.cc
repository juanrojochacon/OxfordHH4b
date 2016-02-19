//samples.cc
#include "samples.h"
#include "run.h"
#include <exception>

using namespace std;

// ************************************ Initialisation ************************************

/*
This routine initialises Pythia8
with all the settings for the shower and underlying event
 */
void InitPythia(  runCard const& rc, sampleCard const& sc, uint32_t const& seed, 
                  Pythia8::Pythia& pythiaRun, double& weight_norm )
{
  // Random seed
  pythiaRun.readString("Random:setSeed = on");
  std::ostringstream o;
  o<<"Random:seed = "<<seed;
  std::cout<<o.str()<<std::endl;
  pythiaRun.readString(o.str());

  // Initialize Les Houches Event File run. List initialization information.
  pythiaRun.readString("Beams:frameType = 4"); 
  
  // Switch off counter
  pythiaRun.readString("Next:numberCount = 0");

  // The shower is QCD only, no QED or weak effects included
  pythiaRun.readString("SpaceShower:QEDshowerByQ  = off"); // QED shower off
  pythiaRun.readString("SpaceShower:QEDshowerByL  = off"); // QED shower off
  pythiaRun.readString("TimeShower:QEDshowerByQ = off");  // QED off on ISR / quarks irradiate photons
  pythiaRun.readString("TimeShower:QEDshowerByL = off");  // QED off on ISR / leptons irradiate photons  
  pythiaRun.readString("TimeShower:QEDshowerByGamma = off");  // Allow photons to branch into lepton or quark pairs 

  // Initial and final state radiation 
  if (rc.pythiaShower)
  {
    pythiaRun.readString("PartonLevel:ISR = on");  // Shower on
    pythiaRun.readString("PartonLevel:FSR = on");  // Shower on
  }
  else
  {
    pythiaRun.readString("PartonLevel:ISR = off"); 
    pythiaRun.readString("PartonLevel:FSR = off");
    pythiaRun.readString("PartonLevel:Remnants = off"); // Disable beam-remnants
  }

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
  std::string ofile;
  ofile="Beams:LHEF = "+sc.eventpath;
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
  weight_norm = pythia_wgt/((double) sc.nevt_sample);      // Unit weight

  return;  
}


void InitHepMC( runCard const& rc, sampleCard const& sc, double& weight_norm)
{
  double gen_xsec = 0;
  double sum_weights = 0;

  std::ifstream hepmc_is( sc.eventpath.c_str() );                   // HepMC input
  for (int iEvent = 0; iEvent < sc.nevt_sample; iEvent++)
    {
      if (! hepmc_is )
      {
        cerr << "Error: HepMC end of file!"<<endl;
        exit(-1);
      }

      HepMC::GenEvent evt;
      evt.read( hepmc_is );
      if( !evt.is_valid() ) 
      {
        cerr << "Error: Invalid HepMC event!" << endl;
        exit(-1);
      }

      if (iEvent % 10000 == 0)
        cout << "Examined "<< iEvent<< " events ... "<<endl;

      // Update to generated xsec
      gen_xsec = evt.cross_section()->cross_section();
      sum_weights += evt.weights()[0];
    } 

  weight_norm = ( gen_xsec / sum_weights );

  // Reset istream
  hepmc_is.close();
}

// ************************************ File Input ************************************

/*
  Get the information on all final state particles
 */
void get_final_state_particles(Pythia8::Pythia& pythiaRun, finalState& particles, double& unit_weight)
{
  unit_weight = 1;

  // Fetch next record
  if (!pythiaRun.next()) 
    if (pythiaRun.info.atEndOfFile()) // Stop showering when the end of the LHE file is reached
    {
      std::cerr << "Info: end of input lhe file reached" << std::endl;
      exit(-1);
    }

  for (int i = 0; i < pythiaRun.event.size(); i++)
  {
    // Get PDG ID
    const int particle_id = pythiaRun.event[i].id();

    // Consider only final state particles
    if( pythiaRun.event[i].status() <= 0 ) continue;
    
    // Get the particle kinematics
    const double E = pythiaRun.event[i].e();
    const double px= pythiaRun.event[i].px();
    const double py= pythiaRun.event[i].py();
    const double pz = pythiaRun.event[i].pz();

    // Form PseudoJet
    fastjet::PseudoJet jet(px,py,pz,E);
    jet.set_user_index(particle_id);

    particles.push_back( jet );
  } // End loop over particles in event
}




/*
  Get the information on all final state particles - HepMC version
 */
void get_final_state_particles(std::ifstream& hepmc_is, finalState& particles, double& unit_weight)
{
  if ( !hepmc_is ) 
  {
    cerr << "Info: end of input hepmc file reached" << endl;
    exit(-1);
  }

  HepMC::GenEvent event;
  event.read( hepmc_is );

  if( !event.is_valid() ) 
  {
    cerr << "Error: Invalid HepMC event!" << endl;
    exit(-1);
  }

  // Unit weight
  unit_weight = event.weights()[0];

  // Conversion factors
  const double momConv = HepMC::Units::conversion_factor(event.momentum_unit(), HepMC::Units::GEV);

  for ( HepMC::GenEvent::particle_iterator p = event.particles_begin();
        p != event.particles_end(); ++p ) 
     if ( !(*p)->end_vertex() && (*p)->status()==1 ) // Is final-state
     {
        HepMC::GenParticle* gp = *p;

        // Particle kinematics
        const double E = momConv*gp->momentum().e();
        const double px = momConv*gp->momentum().px();
        const double py = momConv*gp->momentum().py();
        const double pz = momConv*gp->momentum().pz();

        const int pdg = gp->pdg_id();

        // Form pseudojet
        fastjet::PseudoJet jet(px,py,pz,E);
        jet.set_user_index(pdg);

        // push back
        particles.push_back( jet );
     }
}


