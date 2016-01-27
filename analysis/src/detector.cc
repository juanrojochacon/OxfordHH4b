#include "detector.h"
#include "run.h"
#include "utils.h"
#include "samples.h"

#include <math.h> 

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/SoftKiller.hh"

using namespace std;

// Detector resolution
static const double phiRes=0.1;
static const double etaRes=0.1;

// SoftKiller
static const fastjet::contrib::SoftKiller soft_killer(2.5, 0.4);

// Pythia MinBias
Pythia8::Pythia pythiaPileup(std::string(PYTHIADIR));
bool pythiaInit = false;
void initPythiaPileup()
{
  pythiaPileup.readString("Next:numberShowInfo = 0");
  pythiaPileup.readString("Next:numberShowProcess = 0");
  pythiaPileup.readString("Next:numberShowEvent = 0");
  pythiaPileup.readString("Next:numberCount = 0");

  pythiaPileup.readString("Random:setSeed = on");
  std::ostringstream o;
  o<<"Random:seed = "<<int(pythiaSeed()+34267);
  pythiaPileup.readString(o.str());

  // No hadronization
  pythiaPileup.readString("HadronLevel:all = off"); // Of hadronization

  pythiaPileup.readString("SoftQCD:all = on");
  pythiaPileup.settings.parm("Beams:eCM", 14000);
  pythiaPileup.init();

  pythiaInit = true;
}

void AddPileup( int const& nPileup, finalState& particles )
{
	if (!pythiaInit) initPythiaPileup();
	double dummy;
    for (int iPileup = 0; iPileup < nPileup; ++iPileup)
    	get_final_state_particles(pythiaPileup, particles, dummy);
}

void DetectorSim(finalState input, finalState& output)
{
	AddPileup(npileupEvents(), input);
	if (softKillered()) input = soft_killer(input);
	for(size_t i=0; i<input.size(); i++)
	{
		// Detector granularity
		const double newEta = floor(input[i].eta()/etaRes)*etaRes + etaRes/2.0;
		const double newPhi = floor(input[i].phi()/phiRes)*phiRes + phiRes/2.0;

		// Lengthwise gaussian smear
		const double sm = box_muller(1.0,0.01*GetESmear());

		// Reconstruct smeared jet
		const double pT = sm*input[i].pt();
		const double px = pT*cos(newPhi);
		const double py = pT*sin(newPhi);
		const double pz = pT*sinh(newEta);
		const double E = sqrt(input[i].m2()+px*px + py*py + pz*pz);

 		// Form PseudoJet
    	fastjet::PseudoJet jet(px, py, pz, E);
    	jet.set_user_index(input[i].user_index());

    	output.push_back(jet);
	}
}
