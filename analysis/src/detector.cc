#include "detector.h"
#include "run.h"
#include "utils.h"

#include <math.h> 

#include "fastjet/PseudoJet.hh"

using namespace std;

void DetectorSim(finalState const& input, finalState& output)
{
	const double phiRes=0.1;
	const double etaRes=0.1;

	for(size_t i=0; i<input.size(); i++)
	{
		//cout << "Particle: "<<i <<" eta: "<<input[i].eta()<<" phi: "<<input[i].phi()<<" mass: "<< input[i].m()<<endl;

		const double newEta = floor(input[i].eta()/etaRes)*etaRes + etaRes/2.0;
		const double newPhi = floor(input[i].phi()/phiRes)*phiRes + phiRes/2.0;

		//cout << "New eta: "<<newEta<<" phi: "<<newPhi<<endl;

		const double pT = input[i].pt();
		const double px = pT*cos(newPhi);
		const double py = pT*sin(newPhi);
		const double pz = pT*sinh(newEta);
		const double E = sqrt(input[i].m2()+px*px + py*py + pz*pz);

		// Lengthwise gaussian smear
		const double sm = box_muller(1.0,0.01*GetESmear());

 		// Form PseudoJet
    	fastjet::PseudoJet jet(sm*px, sm*py, sm*pz, sm*E);
    	jet.set_user_index(input[i].user_index());

    	if (input[i].has_user_info())
    		jet.set_user_info(new JetInfo(input[i].user_info<JetInfo>()));
    	
    	output.push_back(jet);

		//cout << "Check: eta: "<<jet.eta()<<" phi: "<<jet.phi()<<" mass: "<< jet.m()<<endl<<endl;
	}


   // 

}