// samples.cc
#include "samples.h"
#include "run.h"

#include "sherpa-xs.hh"
#include <algorithm>
#include <cstdlib>
#include <exception>
#include <iterator>
#include <vector>

using namespace std;

int ceil_div(int x, int y) {
    std::div_t division = std::div(x, y);
    return division.quot + (division.rem ? 1 : 0);
}

// ************************************ Initialisation ************************************

/*
This routine initialises Pythia8
with all the settings for the shower and underlying event
 */
void InitPythia(runCard const& rc, sampleCard const& sc, uint32_t const& seed,
                Pythia8::Pythia& pythiaRun, double& weight_norm) {
    // Random seed
    pythiaRun.readString("Random:setSeed = on");
    std::ostringstream o;
    o << "Random:seed = " << seed;
    std::cout << o.str() << std::endl;
    pythiaRun.readString(o.str());

    // Initialize Les Houches Event File run. List initialization information.
    pythiaRun.readString("Beams:frameType = 4");

    // Switch off counter
    pythiaRun.readString("Next:numberCount = 0");

    // The shower is QCD only, no QED or weak effects included
    pythiaRun.readString("SpaceShower:QEDshowerByQ  = off"); // QED shower off
    pythiaRun.readString("SpaceShower:QEDshowerByL  = off"); // QED shower off
    pythiaRun.readString(
          "TimeShower:QEDshowerByQ = off"); // QED off on ISR / quarks irradiate photons
    pythiaRun.readString(
          "TimeShower:QEDshowerByL = off"); // QED off on ISR / leptons irradiate photons
    pythiaRun.readString("TimeShower:QEDshowerByGamma = off"); // Allow photons to branch into
                                                               // lepton or quark pairs

    // Initial and final state radiation
    if (rc.pythiaShower) {
        pythiaRun.readString("PartonLevel:ISR = on"); // Shower on
        pythiaRun.readString("PartonLevel:FSR = on"); // Shower on
    }
    else {
        pythiaRun.readString("PartonLevel:ISR = off");
        pythiaRun.readString("PartonLevel:FSR = off");
        pythiaRun.readString("PartonLevel:Remnants = off"); // Disable beam-remnants
    }

    // No hadronization
    pythiaRun.readString("HadronLevel:all = off"); // Of hadronization

    // For the time being no  UE or PU included
    pythiaRun.readString("PartonLevel:MPI = off"); // Off multiple interactions (UE)

    // Higgs decays always into bbbar
    // Need to correct by hand the xsecs for the BR(HH->4b) branching fraction
    pythiaRun.readString("25:onMode = off");
    pythiaRun.readString("25:onIfAll = 5 -5");

    // Z always decays into bbbar
    // Need to correct by hand the xsecs for the BR(Z->bb) branching fraction
    pythiaRun.readString("23:onMode = off");
    pythiaRun.readString("23:onIfAll = 5 -5");

    // W always decays hadronically
    // Need to correct by hand the xsecs for the BR(W->qq) branching fraction
    pythiaRun.readString("24:onMode = on");
    pythiaRun.readString("24:offIfAny = 11 12 13 14 15 16");

    // b quarks and do not decay
    // They are treated as stable particles in the detector
    pythiaRun.readString("5:mayDecay = no");
    pythiaRun.readString("-5:mayDecay = no");

    // Read the Les Houches Event File
    std::string ofile;
    ofile = "Beams:LHEF = " + sc.eventpath;
    pythiaRun.readString(ofile.c_str());

    // Main initialization
    pythiaRun.init();

    // Verify pythia event sample
    if (pythiaRun.info.nProcessesLHEF() != 1) {
        std::cerr << "Error: number of subprocesses in LHE sample is not equal to 1 " << std::endl;
        std::cerr << "This code does not support this at the moment" << std::endl;
        exit(-1);
    }

    // Event weight information from Pythia
    const double pythia_wgt = pythiaRun.info.sigmaLHEF(0);           // Total sample weight
    weight_norm             = pythia_wgt / ((double)sc.nevt_sample); // Unit weight

    return;
}

void InitHepMC(runCard const& rc, sampleCard const& sc, double& weight_norm,
	       int& evts_per_subsample, std::vector<long long>& subsample_indices) {
    double comp_xsec = 0;
    double gen_xsec = 0;
    //double weight_norm in args
    int nevts = 0;
    std::cout << "Events are at " << sc.eventpath << "\n";

    std::ifstream hepmc_index(sc.eventpath + ".index");
    if(!hepmc_index) { 
        std::cerr << "Failed to open index file ("<< sc.eventpath << ".index).\n"
		  << "Have you run the indexer (HepMCParser) over your sample?\n";
	std::exit(EXIT_FAILURE);
    }
    hepmc_index >> comp_xsec >> gen_xsec >> weight_norm >> nevts >> evts_per_subsample;

    std::cout << std::scientific
	      << "Computed xsec: " << comp_xsec << "\n"
              << "Generated xsec: " << gen_xsec << "\n"
              << "Unit weight: " << weight_norm << "\n"
	      << std::endl;
    if(sc.nevt_sample != nevts) {
        std::cerr << "nevts in sample card does not match sample index.\n"
		  << "nevts should be " << nevts << "\n";
	std::exit(EXIT_FAILURE);
    }

    if(rc.sub_samplesize % evts_per_subsample != 0) {
	std::cerr << "Subsample size in run card is not a multiple of that in sample index.\n"
		  << "Set subsample size to a multiple of " << evts_per_subsample << ".\n";
	std::exit(EXIT_FAILURE);
    }
   
    std::copy(std::istream_iterator<long long>(hepmc_index),
              std::istream_iterator<long long>(),
	      std::back_inserter(subsample_indices));

    if(subsample_indices.size() != ceil_div(nevts, evts_per_subsample)) {
        std::cerr << "Subsample indices vector has the wrong size: " << subsample_indices.size()
		  << " instead of " << ceil_div(nevts, evts_per_subsample) << "\n"
		  << "Is the index corrupt?\n";
	std::exit(EXIT_FAILURE);
    }
}

// ************************************ File Input ************************************

/*
  Get the information on all final state particles
 */
void get_final_state_particles(Pythia8::Pythia& pythiaRun, finalState& particles,
                               double& unit_weight) {
    unit_weight = 1;

    // Fetch next record
    if (!pythiaRun.next())
        if (pythiaRun.info.atEndOfFile()) // Stop showering when the end of the LHE file is reached
        {
            std::cerr << "Info: end of input lhe file reached" << std::endl;
            exit(-1);
        }

    for (int i = 0; i < pythiaRun.event.size(); i++) {
        // Get PDG ID
        const int particle_id = pythiaRun.event[i].id();

        // Consider only final state particles
        if (pythiaRun.event[i].status() <= 0) continue;

        // Get the particle kinematics
        const double E  = pythiaRun.event[i].e();
        const double px = pythiaRun.event[i].px();
        const double py = pythiaRun.event[i].py();
        const double pz = pythiaRun.event[i].pz();

        // Form PseudoJet
        fastjet::PseudoJet jet(px, py, pz, E);
        jet.set_user_index(particle_id);

        particles.push_back(jet);
    } // End loop over particles in event
}

/*
  Get the information on all final state particles - HepMC version
 */
void get_final_state_particles(std::ifstream& hepmc_is, finalState& particles,
                               double& unit_weight, bool print) {
    if (!hepmc_is) {
        cerr << "Info: end of input hepmc file reached" << endl;
        exit(-1);
    }

    HepMC::GenEvent event;
    event.read(hepmc_is);

    if (!event.is_valid()) {
        cerr << "Error: Invalid HepMC event!" << endl;
        exit(-1);
    }

    if(print) event.print();
    // Unit weight
    unit_weight = event.weights()[0];

    // Conversion factors
    const double momConv =
          HepMC::Units::conversion_factor(event.momentum_unit(), HepMC::Units::GEV);

    for (HepMC::GenEvent::particle_iterator p = event.particles_begin(); p != event.particles_end();
         ++p)
        if (!(*p)->end_vertex() && (*p)->status() == 1) // Is final-state
        {
            HepMC::GenParticle* gp = *p;

            // Particle kinematics
            const double E  = momConv * gp->momentum().e();
            const double px = momConv * gp->momentum().px();
            const double py = momConv * gp->momentum().py();
            const double pz = momConv * gp->momentum().pz();

            const int pdg = gp->pdg_id();

            // Form pseudojet
            fastjet::PseudoJet jet(px, py, pz, E);
            jet.set_user_index(pdg);

            // push back
            particles.push_back(jet);
        }
}
