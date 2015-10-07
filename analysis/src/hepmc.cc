#include "hepmc.h"
#include "run.h"

using namespace std;

void InitHepMC( std::string const& eventfile, int const& nevt_max, double& weight_norm)
{
  double gen_xsec = 0;
  double sum_weights = 0;

  std::ifstream hepmc_is( eventfile.c_str() );            // HepMC input
  for (int iEvent = 0; iEvent < nevt_max; iEvent++)
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

        // Form user_info - should set charge and event ID here
        JetInfo* user_info = new JetInfo(0, 0, 0, 0);
        jet.set_user_info(user_info);

        // push back
        particles.push_back( jet );
     }

  // Verify final state four momenta
  Analysis::VerifyFourMomentum(particles);
}

