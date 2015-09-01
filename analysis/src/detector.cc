#include "detector.h"
#include "run.h"
#include "utils.h"
#include "hepmc.h"

#include <math.h> 

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/SoftKiller.hh"

using namespace std;

static std::ifstream pileupStream;
static int pileupCount = 0;

// SoftKiller
static fastjet::contrib::SoftKiller soft_killer(3, 0.1);

void AddPileup( int const& nPileup, finalState& particles )
{
	for ( int iEvent = 0; iEvent < nPileup; iEvent++ )
	{
		 if (!pileupStream || pileupCount >= npileupTotal())
		 {
		 	pileupStream.clear() ;
			pileupStream.seekg(0, ios::beg);
			pileupCount = 0;
		 }

		 if (!pileupStream.is_open())
		 	pileupStream.open( std::string(SAMPLEDIR) + minBiasFile() );

		 double dummy;
		 get_final_state_particles(pileupStream, particles, dummy);
		 pileupCount++;
	}
}

typedef std::map<std::pair<int, int>, fastjet::PseudoJet > JetMap;
void DetectorSim(finalState input, finalState& output)
{
	if (pileupSimulated())
	{
		AddPileup(npileupEvents(), input);
		input = soft_killer(input);
	}

	const double phiRes=0.1;
	const double etaRes=0.1;

	//Form Hadronic calorimeter towers
	JetMap caloTowers;
	for(size_t i=0; i<input.size(); i++)
	{
		const double newEta = floor(input[i].eta()/etaRes)*etaRes + etaRes/2.0;
		const double newPhi = floor(input[i].phi()/phiRes)*phiRes + phiRes/2.0;

		const double pT = input[i].pt();
		const double px = pT*cos(newPhi);
		const double py = pT*sin(newPhi);
		const double pz = pT*sinh(newEta);
		const double E = sqrt(input[i].m2()+px*px + py*py + pz*pz);

 		// Form PseudoJet
    	fastjet::PseudoJet jet(px, py, pz, E);
    	jet.set_user_index(input[i].user_index());
    	if (input[i].has_user_info())
    		jet.set_user_info(new JetInfo(input[i].user_info<JetInfo>()));

    	// Look for existing calo tower
    	const std::pair<int,int> caloCoords(newEta*100.0, newPhi*100.0);
    	JetMap::iterator f = caloTowers.find(caloCoords);
    	if (f != caloTowers.end())
    	{
    		(*f).second += jet;

    		// Set highest user_index
    		const int jidx = jet.user_index();
    		if (abs(jidx) < 6 && jet.pt() > 15)
    			if (abs(jidx) > abs((*f).second.user_index()))
    				(*f).second.set_user_index(jidx);

    	}
    	else
    	{
    		caloTowers[caloCoords] = jet;
    	}
	}

	for(JetMap::iterator iTower = caloTowers.begin(); iTower != caloTowers.end(); iTower++)
	{
		const fastjet::PseudoJet& tower = (*iTower).second;

		// Lengthwise gaussian smear
		const double sm = box_muller(1.0,0.01*GetESmear());
		const double px = sm*tower.px();
		const double py = sm*tower.py();
		const double pz = sm*tower.pz();

		// Tower is massless
		const double E = sqrt(px*px + py*py + pz*pz);

 		// Form PseudoJet
    	fastjet::PseudoJet jet(px, py, pz, E);
    	jet.set_user_index((*iTower).second.user_index());
    	output.push_back(jet);
	}
}