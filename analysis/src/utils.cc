#include "utils.h"

#include <cmath>
#include <cstdlib>
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

using namespace std;

double getDPhi(double phi1, double phi2) {
    double deltaPhi                = phi1 - phi2;
    if (deltaPhi > M_PI) deltaPhi  = (deltaPhi - 2.0 * M_PI);
    if (deltaPhi < -M_PI) deltaPhi = (deltaPhi + 2.0 * M_PI);
    return deltaPhi;
}

// ----------------------------------------------------------------------------------
// Recluster with kt algorithm to obtain splitting scales
double SplittingScales(fastjet::PseudoJet const& jet) {
    if (!jet.has_constituents()) {
        std::cerr << "ERROR! Splittings (d12, d23, ...) can only be calculated on jets for which "
                     "the constituents are known."
                  << std::endl;
        return -1;
    }

    vector<fastjet::PseudoJet> const& constits = jet.constituents();

    if (constits.size() == 0) {
        std::cerr << "ERROR! Splittingscales must have >0 constituents!" << std::endl;
        return -1;
    }

    fastjet::JetDefinition ekt_jd =
          fastjet::JetDefinition(fastjet::kt_algorithm, 1.5, fastjet::E_scheme, fastjet::Best);
    const fastjet::ClusterSequence  kt_seq_excl  = fastjet::ClusterSequence(constits, ekt_jd);
    std::vector<fastjet::PseudoJet> excl_cluster = sorted_by_pt(kt_seq_excl.inclusive_jets());

    if (excl_cluster.size() == 0) {
        std::cerr << "ERROR! Splittingscales must have >0 excl_cluster!" << std::endl;
        return -1;
    }

    fastjet::PseudoJet kt_jet = excl_cluster[0];

    const double split12 = 1.5 * sqrt(kt_seq_excl.exclusive_subdmerge(kt_jet, 1));
    return split12;
}

std::vector<double> SplittingScales(std::vector<fastjet::PseudoJet> const& jetVec) {
    // vectors that contain the respective splitting scales for all jets
    std::vector<double> split12_vec;

    for (const auto& i : jetVec) split12_vec.push_back(SplittingScales(i));

    return split12_vec;
}

// ----------------------------------------------------------------------------------
// Calculate jet pull vector
std::vector<double> JetPullVector(fastjet::PseudoJet const& jet) {

    std::vector<double> jetPullVector;
    jetPullVector.push_back(-1);
    jetPullVector.push_back(-1);

    if (!jet.has_constituents()) {
        std::cerr << "ERROR! Jet pull vector can only be calculated on jets for which the "
                     "constituents are known."
                  << std::endl;
        return jetPullVector;
    }

    std::vector<fastjet::PseudoJet> const& constituents = jet.constituents();

    if (constituents.size() == 0) {
        std::cerr << "ERROR! Jet pull vector can only be calculated from jets with >0 constituents!"
                  << std::endl;
        return jetPullVector;
    }

    double jetPull_y   = 0;
    double jetPull_phi = 0;

    // std::cout << "Number of constituents " << constituents.size() << std::endl;

    for (const auto& constit : constituents) {

        double dY      = constit.rapidity() - jet.rapidity();
        double dPhi    = getDPhi(constit.phi(), jet.phi());
        double absDiff = sqrt(pow(dY, 2) + pow(dPhi, 2));

        double pTRatio = constit.pt() * pow(jet.pt(), -1);
        // std::cout << "dY " << dY << std::endl;
        // std::cout << "dPhi " << dPhi << std::endl;
        // std::cout << "absDiff " << absDiff << std::endl;
        // std::cout << "constit.pt() " << constit.pt() << std::endl;
        // std::cout << "jet.pt() " << jet.pt() << std::endl;
        // std::cout << "pTRatio " << pTRatio << std::endl;

        jetPull_y += pTRatio * absDiff * dY;
        jetPull_phi += pTRatio * absDiff * dPhi;
    }

    while (jetPull_phi > M_PI) jetPull_phi  = (jetPull_phi - 2.0 * M_PI);
    while (jetPull_phi < -M_PI) jetPull_phi = (jetPull_phi + 2.0 * M_PI);

    jetPullVector[0] = jetPull_y;
    jetPullVector[1] = jetPull_phi;

    // std::cout << "jetPull_y " << jetPull_y << std::endl;
    // std::cout << "jetPull_phi " << jetPull_phi << std::endl;

    return jetPullVector;
}

// ----------------------------------------------------------------------------------
// Calculate jet pull of jet1 w.r.t jet2
double JetPull(fastjet::PseudoJet const& jet1, fastjet::PseudoJet const& jet2) {

    std::vector<double> jetPullVector = JetPullVector(jet1);

    if (jetPullVector[0] == -1 || jetPullVector[1] == -1) {
        std::cerr << "ERROR! Jet pull vector for jet1 cannot be calculated." << std::endl;
        return -1;
    }

    fastjet::PseudoJet jetSep = jet1 - jet2;

    double lenPullVec = jetPullVector[0] * jetPullVector[0] + jetPullVector[1] * jetPullVector[1];
    double lenjetSep  = jetSep.rapidity() * jetSep.rapidity() + jetSep.phi() * jetSep.phi();

    if (lenPullVec == 0.) {
        std::cerr << "ERROR! Jet pull vector has zero length." << std::endl;
        return -1;
    }

    if (lenjetSep == 0.) {
        std::cerr << "ERROR! Jet separation vector has zero length." << std::endl;
        return -1;
    }

    double dotProd = jetPullVector[0] * jetSep.rapidity() + jetPullVector[1] * jetSep.phi();
    dotProd /= sqrt(lenPullVec);
    dotProd /= sqrt(lenjetSep);

    double jetPull = std::acos(dotProd);

    // std::cout << "dotProd " << dotProd << std::endl;
    // std::cout << "jetPull " << jetPull << std::endl;

    return jetPull;
}

// ----------------------------------------------------------------------------------
// Calculate angular variables
double Chi(const fastjet::PseudoJet& h0, const fastjet::PseudoJet& h1) {

    double y0 = h0.rapidity();
    double y1 = h1.rapidity();
    return exp(fabs(y0 - y1));
}

// ---------------------------------------------------------------------------------------
// Associate charged small-R ("track") jets to large-R ("calo") jet
void get_assoc_trkjets(const fastjet::PseudoJet& calojet, std::vector<fastjet::PseudoJet> trkjets,
                       std::vector<fastjet::PseudoJet>& matched_trkjets, bool debug = false) {

    // vector to hold input clusters and ghosts
    std::vector<fastjet::PseudoJet> input_particles;
    input_particles.clear();

    // jet clusters from large-R jet
    vector<fastjet::PseudoJet> constituents = calojet.constituents();

    if (debug) std::cout << "calo constituents size = " << constituents.size() << std::endl;

    for (auto noghost : constituents) {

        noghost.reset_PtYPhiM(noghost.pt(), noghost.rapidity(), noghost.phi(), 0.0);
        if (noghost.E() < 0.) continue;

        // set user index for calo clusters to -1 to differentiate them from "track" constituents
        // later
        noghost.set_user_index(-1);
        input_particles.push_back(noghost);
    }

    if (debug)
        std::cout << "calo only input particles size = " << input_particles.size() << std::endl;

    // make ghost PseudoJets out of track jet direction
    for (unsigned int trackJetItr = 0; trackJetItr < trkjets.size(); ++trackJetItr) {

        fastjet::PseudoJet myghost = trkjets.at(trackJetItr);
        // if( myghost.pt() <= 20.0 || fabs( myghost.rapidity() ) >= 2.5 ) continue;

        myghost.reset_PtYPhiM(1e-12, myghost.rapidity(), myghost.phi(), 0.0);
        if (myghost.E() < 0.) continue;

        myghost.set_user_index(trackJetItr);
        input_particles.push_back(myghost);
    }

    if (debug)
        std::cout << "calo+track jets input particles size = " << input_particles.size()
                  << std::endl;

    // do ghost association and get list of pseudojet track jets that are associated
    double                       Rparam        = 1.0;
    fastjet::Strategy            strategy      = fastjet::Best;     // according to atlas reco
    fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme; // according to atlas reco
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, recomb_scheme, strategy);

    // run the jet clustering with the above jet definition
    fastjet::ClusterSequence   clust_seq(input_particles, jet_def);
    vector<fastjet::PseudoJet> sorted_jets = fastjet::sorted_by_pt(clust_seq.inclusive_jets());

    if (debug) std::cout << "number of sorted jets = " << sorted_jets.size() << std::endl;

    fastjet::PseudoJet newJet =
          sorted_jets.at(0); // there are more jets in the vector, but they all have pT ~0
    if (debug)
        std::cout << "new jet constituent size = " << newJet.constituents().size() << std::endl;
    vector<fastjet::PseudoJet> newJet_constituents = newJet.constituents();
    for (const auto& constit : newJet_constituents) {
        if (debug)
            std::cout << " user index = " << constit.user_index()
                      << ", pt of constit = " << constit.pt() << std::endl;
        if (constit.user_index() >= 0) {
            int iter = constit.user_index();
            // 	  if(trkjets.at(iter).pt() > 20. && fabs(trkjets.at(iter).eta()) < 2.5 )
            // matched_trkjets.push_back(trkjets.at(iter));
            matched_trkjets.push_back(trkjets.at(iter));
        }
    }

    // Sort matched jets by pt
    matched_trkjets = sorted_by_pt(matched_trkjets);

    return;
}

// ---------------------------------------------------------------------------------------
// B-tagging based on ATLAS PubNotes
// (ATL-PHYS-PUB-2014-013, ATL-PHYS-PUB-2014-014)

double btag_eff(double jet_pt) {

    if (jet_pt <= 20.)
        return 0.05;
    else if (jet_pt > 20. && jet_pt <= 50.)
        return 0.35;
    else if (jet_pt > 50. && jet_pt <= 70.)
        return 0.65;
    else if (jet_pt > 70. && jet_pt <= 100.)
        return 0.70;
    else if (jet_pt > 100. && jet_pt <= 200.)
        return 0.65;
    else if (jet_pt > 200. && jet_pt <= 350.)
        return 0.60;
    else if (jet_pt > 350. && jet_pt <= 500.)
        return 0.55;
    else
        return 0.50;
}

double mistag_eff(double jet_pt) {

    if (jet_pt <= 70.)
        return 1. / 150.;
    else if (jet_pt > 70. && jet_pt <= 180.)
        return 1. / 170.;
    else if (jet_pt > 180. && jet_pt <= 250.)
        return 1. / 130.;
    else
        return 1. / 100.;
}

double charm_eff(double jet_pt) {

    if (jet_pt <= 70.)
        return 1. / 6.;
    else if (jet_pt > 70. && jet_pt <= 180.)
        return 1. / 4.;
    else if (jet_pt > 180. && jet_pt <= 250.)
        return 1. / 4.5;
    else
        return 1. / 5.;
}
