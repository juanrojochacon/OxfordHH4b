#include "detector.h"
#include "run.h"
#include "utils.h"

#include <math.h> 

#include "fastjet/PseudoJet.hh"

using namespace std;

typedef std::map<std::pair<int, int>, fastjet::PseudoJet > JetMap;
void DetectorSim(finalState const& input, finalState& output)
{
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