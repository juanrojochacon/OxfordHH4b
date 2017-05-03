#include "detector.h"
#include "run.h"
#include "samples.h"
#include "utils.h"

#include <cmath>
#include <random>

#include "fastjet/PseudoJet.hh"

using namespace std;

// Detector simulation

Detector::Detector(runCard const& run, sampleCard const& sample, int const& pythiaSeed,
                   int const& detectorSeed)
    : pythiaPileup(std::string(PYTHIADIR)),
      rng(detectorSeed),
      nPileup(run.npileup),
      phiRes(0.1),
      etaRes(0.1),
      jetEsmear(run.jetEsmear) {
    pythiaPileup.readString("Next:numberShowInfo = 0");
    pythiaPileup.readString("Next:numberShowProcess = 0");
    pythiaPileup.readString("Next:numberShowEvent = 0");
    pythiaPileup.readString("Next:numberCount = 0");

    pythiaPileup.readString("Random:setSeed = on");
    std::ostringstream o;
    o << "Random:seed = " << pythiaSeed;
    pythiaPileup.readString(o.str());

    // No hadronization
    pythiaPileup.readString("HadronLevel:all = off"); // Of hadronization

    pythiaPileup.readString("SoftQCD:all = on");
    pythiaPileup.settings.parm("Beams:eCM", sample.sqrts);
    pythiaPileup.init();
}

void Detector::AddPileup(finalState& particles) {
    double dummy;
    for (int iPileup = 0; iPileup < nPileup; ++iPileup)
        get_final_state_particles(pythiaPileup, particles, dummy);
}

void Detector::Simulate(finalState input, finalState& output) {
    AddPileup(input);
    for (auto& i : input) {
        // Detector granularity
        const double newEta = floor(i.eta() / etaRes) * etaRes + etaRes / 2.0;
        const double newPhi = floor(i.phi() / phiRes) * phiRes + phiRes / 2.0;

        // Lengthwise gaussian smear
        std::normal_distribution<> normal_dist(1.0, 0.01 * jetEsmear);
        const double               sm = normal_dist(rng);

        // Reconstruct smeared jet
        const double pT = sm * i.pt();
        const double px = pT * cos(newPhi);
        const double py = pT * sin(newPhi);
        const double pz = pT * sinh(newEta);
        const double E  = sqrt(i.m2() + px * px + py * py + pz * pz);

        // Form PseudoJet
        fastjet::PseudoJet jet(px, py, pz, E);
        jet.set_user_index(i.user_index());

        output.push_back(jet);
    }
}

/*

bool Analysis::VerifyFourMomentum(std::vector<fastjet::PseudoJet> const& jets)
{
  // Smearing breaks four-mom verification
  if ( GetESmear() > 1E-8 ) return true;
  // No beam-remnants
  if (!pythiaShowered()) return true;


  // To check energy-momentum conservation
  double const Eref=14000; // Samples generated for LHC 14 TeV
  double const tol_emom=1.0;

  // Check again four-momentum conservation, this time applied to jets
  // formed from the clustering of quarks and gluons (and beam remnants as well)
  double px_tot=0;
  double py_tot=0;
  double pz_tot=0;
  double E_tot=0;
  for(size_t ij=0;ij<jets.size();ij++)
  {
    px_tot+= jets[ij].px();
    py_tot+= jets[ij].py();
    pz_tot+= jets[ij].pz();
    E_tot+= jets[ij].E();
  }

  // Check energy-momentum conservation
  if( fabs(px_tot) > tol_emom || fabs(py_tot)  > tol_emom
  || fabs(pz_tot)  > tol_emom || fabs(E_tot-Eref)  > tol_emom )
  {
    std::cout<<"\n **********************************************************************
\n"<<std::endl;
    std::cout<<"No conservation of energy after jet reconstruction "<<std::endl;
    std::cout<<"N_particles = " << jets.size()<<std::endl;
    std::cout<<"px_tot = "<<px_tot<<std::endl;
    std::cout<<"py_tot = "<<py_tot<<std::endl;
    std::cout<<"pz_tot = "<<pz_tot<<std::endl;
    std::cout<<"E_tot, Eref = "<<E_tot<<" "<<Eref<<std::endl;
    std::cout<<"\n **********************************************************************
\n"<<std::endl;
    return false;
  }
  return true;
}
*/
