// oxford_res_fr.cc

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

#include "oxford_sideband.h"
#include "run.h"
#include "utils.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/MassDropTagger.hh"

#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftKiller.hh"

#include <algorithm>

using std::vector;
using namespace fastjet::contrib;
using namespace fastjet;

// LST energy correlations
const EnergyCorrelatorC2 C2(2, EnergyCorrelator::pt_R);
const EnergyCorrelatorD2 D2(2, EnergyCorrelator::pt_R);

// tau2/tau1 NSubjettiness
const NsubjettinessRatio tau21(2, 1, KT_Axes(), UnnormalizedMeasure(1));

// Debugging
const bool debug = false;

// Exclusivity cut
const bool exclusive = false;

// Analysis settings
const int         nAnalysis          = 3;
const int         nCuts              = 5;
const std::string aString[nAnalysis] = {"_res", "_inter", "_boost"};
const std::string cString[nCuts]     = {"_GEN", "_RCO", "_SIG", "_SDBA", "_SDBB"};

// **************************** Reconstruction helper functions ****************************

static vector<PseudoJet> trimJets(vector<PseudoJet> const& injets) {
    const double          Rfilt           = 0.2;
    const double          pt_fraction_min = 0.05;
    const fastjet::Filter trimmer(Rfilt, fastjet::SelectorPtFractionMin(pt_fraction_min));
    vector<PseudoJet>     trimmedJets;
    for (size_t i = 0; i < injets.size(); i++) trimmedJets.push_back(trimmer(injets[i]));
    return trimmedJets;
}

// Check if jets are mass-drop tagged
static vector<PseudoJet> MDtagJets(vector<PseudoJet> const& injets) {
    // Mass-drop tagger
    double const mu   = 0.67;
    double const ycut = 0.09;

    const fastjet::JetDefinition  CA10(fastjet::cambridge_algorithm, 1.0);
    const fastjet::MassDropTagger md_tagger(mu, ycut);
    vector<PseudoJet>             MDTJets;
    for (size_t i = 0; i < injets.size(); i++) {
        const fastjet::ClusterSequence cs_sub(injets[i].constituents(), CA10);
        vector<PseudoJet>              ca_jets    = sorted_by_pt(cs_sub.inclusive_jets());
        const PseudoJet                ca_jet     = ca_jets[0];
        const PseudoJet                tagged_jet = md_tagger(ca_jet);
        if (tagged_jet != 0) MDTJets.push_back(injets[i]);
    }
    return MDTJets;
}

// Returns the two hardest GA subjets for each input largeR jet
static vector<vector<PseudoJet>> getSubJets(vector<PseudoJet> const& largeRJets,
                                            vector<PseudoJet> const& trackJets) {
    vector<vector<PseudoJet>> largeRsubJets;
    for (PseudoJet jet : largeRJets) {
        vector<PseudoJet> subJets;
        get_assoc_trkjets(jet, trackJets, subJets, false);
        largeRsubJets.push_back(SelectorNHardest(2)(sorted_by_pt(subJets)));
    }
    return largeRsubJets;
}

// Small-R B-tagging
static vector<btagType> BTagging(vector<PseudoJet> const& jets_vec) {
    vector<btagType> btag_vec;
    for (auto jet : jets_vec) {
        btagType                 type             = NTAG;
        const vector<PseudoJet>& jet_constituents = SelectorPtMin(15)(jet.constituents());
        for (auto constituent : jet_constituents) {
            const int userid                           = constituent.user_index();
            if (abs(userid) == 5) type                 = BTAG;
            if (abs(userid) == 4 && type != BTAG) type = CTAG;
            if (type == NTAG) type                     = LTAG;
        }
        btag_vec.push_back(type);
    }
    return btag_vec;
}

// nTag: How many b-tags are required
// nB: How many true b-jets are present
// nC: How many true c-jets are present
// nL: How many light jets are present
static double btagProb(int const& nTag, int const& nB, int const& nC, int const& nL) {
    // Choose working point with high purity
    const double btag_prob   = 0.80; // Probability of correct b tagging
    const double btag_mistag = 0.01; // Mistag probability
    const double ctag_prob   = 0.1;  // 0.17; // c-mistag rate ~1/6.

    // Probability of all permutations
    double totalProb = 0;

    // Loop over all possible classification permutations
    for (int iB = 0; iB <= std::min(nB, nTag); iB++)                   // iB b-jets tagged as b
        for (int iC = 0; iC <= std::min(nC, nTag - iB); iC++)          // iC c-jets tagged as b
            for (int iL = 0; iL <= std::min(nL, nTag - iB - iC); iL++) // iL l-jets tagged as b
            {
                const double bProb = pow(btag_prob, iB) * pow(1.0 - btag_prob, nB - iB);
                const double cProb = pow(ctag_prob, iC) * pow(1.0 - ctag_prob, nC - iC);
                const double lProb = pow(btag_mistag, iL) * pow(1.0 - btag_mistag, nL - iL);

                // Does the current permutation have the correct number of b-Tags?
                const int permutationTags = iB + iC + iL; // Number of b-Tags in current permutation
                if (permutationTags == nTag) totalProb += bProb * cProb * lProb;
            }

    return totalProb;
}

OxfordSidebandAnalysis::OxfordSidebandAnalysis(runCard const& run, sampleCard const& sample,
                                               int const& subsample)
    : Analysis("sideband", run, sample, subsample), subtractPU(run.npileup > 0) {
    // ********************* Histogram settings******************

    const double DeltaRmin = 0;
    const double DeltaRmax = 5;

    const double DeltaPhimin = -3.2;
    const double DeltaPhimax = 3.2;

    const double DeltaEtamin = -2.5;
    const double DeltaEtamax = 2.5;

    const double m_min = 0.;
    const double m_max = 180.;

    const double pt_min = 0.;
    const double pt_max = 900.;

    const double eta_min = -6.;
    const double eta_max = +6.;

    const double phi_min = -3.15;
    const double phi_max = +3.15;

    const double m_HH_min = 0.;
    const double m_HH_max = 600.;

    const double pt_HH_min = 0.;
    const double pt_HH_max = 300.;

    const double chi_HH_min = 0.;
    const double chi_HH_max = 30.;

    const int nbins = 30;

    // ********************* Histogram definitions ******************

    for (int i = 0; i < nAnalysis; i++) {
        BookHistogram(new YODA::Histo1D(nCuts, 0, nCuts), "CF" + aString[i]);
        BookHistogram(new YODA::Histo1D(nCuts, 0, nCuts), "CFN" + aString[i]);

        for (int j = 0; j < nCuts; j++) {
            const std::string suffix = aString[i] + cString[j];

            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_smallR" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_smallR" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_smallR" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_largeR" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_largeR" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_largeR" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_H0" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_H1" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, phi_min, phi_max), "phi_H0" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, phi_min, phi_max), "phi_H1" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, pt_HH_min, pt_HH_max), "pt_HH" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, m_HH_min, m_HH_max), "m_HH" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, DeltaRmin, DeltaRmax), "dR_HH" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, DeltaPhimin, DeltaPhimax), "dPhi_HH" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, DeltaEtamin, DeltaEtamax), "dEta_HH" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, chi_HH_min, chi_HH_max), "chi_HH" + suffix);

            BookHistogram(new YODA::Histo2D(nbins, pt_min, pt_max, nbins, pt_min, pt_max),
                          "ptHptH" + suffix);
            BookHistogram(new YODA::Histo2D(nbins, m_min, m_max, nbins, m_min, m_max),
                          "mHmH" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_leadSJ_fj1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_subleadSJ_fj1" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_leadSJ_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_subleadSJ_fj2" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_leadSJ_fj" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_subleadSJ_fj" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_leadSJ_fj1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_subleadSJ_fj1" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_leadSJ_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_subleadSJ_fj2" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_leadSJ_fj" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_subleadSJ_fj" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj" + suffix);

            BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj" + suffix);

            BookHistogram(new YODA::Histo1D(nbins * 2, 0, 1), "C2_fj1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins * 2, 0, 1), "C2_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins * 2, 0, 1), "C2_fj" + suffix);

            BookHistogram(new YODA::Histo1D(nbins * 2, 0, 1), "D2_fj1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins * 2, 0, 1), "D2_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins * 2, 0, 1), "D2_fj" + suffix);

            // Additional histograms from unreweighted analysis
            BookHistogram(new YODA::Histo1D(10, 0, 10), "N_SmallRJets" + suffix);
            BookHistogram(new YODA::Histo1D(10, 0, 10), "N_SmallRJets_BJets" + suffix);
            BookHistogram(new YODA::Histo1D(10, 0, 10), "N_SmallRJets_BTagged" + suffix);
        }
    }

    // ********************* Ntuple definition **********************

    const std::string tupleSpec = "# signal source weight pt_H0 pt_H1 pt_HH m_H0 m_H1 m_HH dR_HH "
                                  "dPhi_HH dEta_HH chi_HH pt_H0_sub0 pt_H0_sub1 pt_H1_sub0 "
                                  "pt_H1_sub1";

    const std::string root = "." + GetRoot() + GetSample() + "/";
    std::stringstream suffix;
    suffix << "." << GetSubSample() << ".dat";

    const std::string resDir = root + "resNTuple" + suffix.str();
    const std::string intDir = root + "intNTuple" + suffix.str();
    const std::string bstDir = root + "bstNTuple" + suffix.str();

    resNTuple.open(resDir.c_str());
    intNTuple.open(intDir.c_str());
    bstNTuple.open(bstDir.c_str());

    resNTuple << tupleSpec << std::endl;
    intNTuple << tupleSpec << " split12_fj tau21_fj C2_fj D2_fj" << std::endl;
    bstNTuple << tupleSpec
              << " split12_fj1 split12_fj2 tau21_fj1 tau21_fj2 C2_fj1 C2_fj2 D2_fj1 D2_fj2"
              << std::endl;

    std::cout << "Oxford PU subtraction: " << subtractPU << std::endl;
}

void OxfordSidebandAnalysis::Analyse(bool const& signal, double const& weightnorm,
                                     finalState const& ifs) {
    Analysis::Analyse(signal, weightnorm, ifs);

    // Perform softKiller subtraction
    const fastjet::contrib::SoftKiller soft_killer(2.5, 0.4);
    const finalState                   fs = subtractPU ? soft_killer(ifs) : ifs;

    // *************************************** General cut selectors
    // *************************************

    const Selector LR_kinematics =
          SelectorNHardest(2) * (SelectorAbsRapMax(2.0) && SelectorPtMin(200.0));
    const Selector SR_kinematics =
          SelectorNHardest(4) * (SelectorAbsRapMax(2.5) && SelectorPtMin(40.0));
    const Selector TR_kinematics = SelectorAbsRapMax(2.5) && SelectorPtMin(50.0);

    // ********************************* Jet clustering  ***************************************

    // Cluster small-R jets
    const double                   ResJetR = 0.4;
    const fastjet::JetDefinition   akt_res(fastjet::antikt_algorithm, ResJetR);
    const fastjet::ClusterSequence cs_akt_res(fs, akt_res);
    const vector<PseudoJet> smallRJets = sorted_by_pt(SR_kinematics(cs_akt_res.inclusive_jets()));
    const vector<btagType>  tagType_SR = BTagging(smallRJets);

    // Cluster small-R track jets
    const double                   GAjetR = 0.3; // Boosted subjet radius for ghost-association
    const fastjet::JetDefinition   jd_subjets(fastjet::antikt_algorithm, GAjetR);
    const fastjet::ClusterSequence cs_subjets(fs, jd_subjets);
    const vector<PseudoJet> trackJets = sorted_by_pt(TR_kinematics(cs_subjets.inclusive_jets()));

    // Cluster large-R jets
    const double                   BoostJetR = 1.0;
    const fastjet::JetDefinition   akt_boost(fastjet::antikt_algorithm, BoostJetR);
    const fastjet::ClusterSequence cs_akt_bst(fs, akt_boost);
    const vector<PseudoJet>        largeRJets_noTrim = LR_kinematics(cs_akt_bst.inclusive_jets());
    const vector<PseudoJet>        largeRJets_Trim =
          subtractPU ? trimJets(largeRJets_noTrim) : largeRJets_noTrim;
    const vector<PseudoJet>         largeRJets    = sorted_by_pt(MDtagJets(largeRJets_Trim));
    const vector<vector<PseudoJet>> largeRsubJets = getSubJets(largeRJets, trackJets);
    vector<vector<btagType>>        tagType_LR;
    for (auto subjets : largeRsubJets) tagType_LR.push_back(BTagging(subjets));

    // ***************************************** Initial histograms
    // **********************************************

    FillHistogram("CF_res", weightnorm, 0.1);
    FillHistogram("CF_inter", weightnorm, 0.1);
    FillHistogram("CF_boost", weightnorm, 0.1);

    FillHistogram("CFN_res", 1., 0.1);
    FillHistogram("CFN_inter", 1., 0.1);
    FillHistogram("CFN_boost", 1., 0.1);

    // **************************************** Boosted analysis
    // *********************************************

    // Set initial weight
    if (exclusive) {
        double event_weight = weightnorm;
        event_weight -=
              BoostedAnalysis(largeRJets, largeRsubJets, tagType_LR, signal, event_weight);
        event_weight -= IntermediateAnalysis(largeRJets, smallRJets, largeRsubJets, tagType_LR,
                                             tagType_SR, signal, event_weight);
        event_weight -= ResolvedAnalysis(smallRJets, tagType_SR, signal, event_weight);
    }
    else {
        BoostedAnalysis(largeRJets, largeRsubJets, tagType_LR, signal, weightnorm);
        ResolvedAnalysis(smallRJets, tagType_SR, signal, weightnorm);
        IntermediateAnalysis(largeRJets, smallRJets, largeRsubJets, tagType_LR, tagType_SR, signal,
                             weightnorm);
    }
    return;
}

double OxfordSidebandAnalysis::BoostedAnalysis(vector<PseudoJet> const&         largeRJets,
                                               vector<vector<PseudoJet>> const& largeRsubJets,
                                               vector<vector<btagType>> btags, bool const& signal,
                                               double const& event_weight) {
    if (largeRJets.size() == 2 && btags.size() == largeRJets.size()) {
        // b-Tagging rate
        int nB = 0, nC = 0, nL = 0;
        for (auto tagvec : btags) {
            // These vectors are restricted to the two hardest subjets already
            nB += std::count(tagvec.begin(), tagvec.end(), BTAG);
            nC += std::count(tagvec.begin(), tagvec.end(), CTAG);
            nL += std::count(tagvec.begin(), tagvec.end(), LTAG);
        }

        if (nB + nC + nL == 4) // Selected 4 candidates
        {
            const double selEff = btagProb(4, nB, nC, nL);
            const double selWgt = selEff * event_weight;
            HiggsFill(largeRJets[0], largeRJets[1], "boost", 1, selWgt);

            const fastjet::PseudoJet dihiggs_boost = largeRJets[0] + largeRJets[1];

            // Higgs mass-window
            const double diffHiggs_0 = fabs(largeRJets[0].m() - 125.);
            const double diffHiggs_1 = fabs(largeRJets[1].m() - 125.);
            const double massWindow  = 20.0;

            const bool signal0 = diffHiggs_0 < massWindow;
            const bool signal1 = diffHiggs_1 < massWindow;

            const bool control0 = !signal0 && diffHiggs_0 < 2.0 * massWindow;
            const bool control1 = !signal1 && diffHiggs_1 < 2.0 * massWindow;

            if (signal0 && signal1) {

                HiggsFill(largeRJets[0], largeRJets[1], "boost", 2, selWgt);
                BoostFill(largeRJets[0], largeRJets[1], "boost", 2, selWgt);

                // Calculate some substructure variables
                const std::vector<double> split12_vec = SplittingScales(largeRJets);
                const double              tau21_fj1   = tau21(largeRJets[0]);
                const double              tau21_fj2   = tau21(largeRJets[1]);

                // C2 energy correlation double-ratio
                const double C2_fj1 = C2(largeRJets[0]);
                const double C2_fj2 = C2(largeRJets[1]);

                // D2 energy correlation double-ratio
                const double D2_fj1 = D2(largeRJets[0]);
                const double D2_fj2 = D2(largeRJets[1]);

                // Fill tuple
                bstNTuple << signal << "\t" << GetSample() << "\t" << selWgt << "\t"
                          << largeRJets[0].pt() << "\t" << largeRJets[1].pt() << "\t"
                          << dihiggs_boost.pt() << "\t" << largeRJets[0].m() << "\t"
                          << largeRJets[1].m() << "\t" << dihiggs_boost.m() << "\t"
                          << largeRJets[0].delta_R(largeRJets[1]) << "\t"
                          << getDPhi(largeRJets[0].phi(), largeRJets[1].phi()) << "\t"
                          << fabs(largeRJets[0].eta() - largeRJets[1].eta()) << "\t"
                          << Chi(largeRJets[0], largeRJets[1]) << "\t" << largeRsubJets[0][0].pt()
                          << "\t" << largeRsubJets[0][1].pt() << "\t" << largeRsubJets[1][0].pt()
                          << "\t" << largeRsubJets[1][1].pt() << "\t" << split12_vec[0] << "\t"
                          << split12_vec[1] << "\t" << tau21_fj1 << "\t" << tau21_fj2 << "\t"
                          << C2_fj1 << "\t" << C2_fj2 << "\t" << D2_fj1 << "\t" << D2_fj2 << "\t"
                          << std::endl;
            }
            else if ((signal0 && control1) || (signal1 && control0)) {
                HiggsFill(largeRJets[0], largeRJets[1], "boost", 3, selWgt);
                BoostFill(largeRJets[0], largeRJets[1], "boost", 3, selWgt);
            }
            else if (control0 && control1) {
                HiggsFill(largeRJets[0], largeRJets[1], "boost", 4, selWgt);
                BoostFill(largeRJets[0], largeRJets[1], "boost", 4, selWgt);
            }
            return selWgt;
        }
    }
    return 0;
};

double OxfordSidebandAnalysis::ResolvedAnalysis(vector<PseudoJet> const& srj,
                                                vector<btagType> const& btags, bool const& signal,
                                                double const& event_weight) {
    if (srj.size() >= 4) {
        // Reconstruct Higgs candidates from small-R jets
        typedef std::pair<PseudoJet, PseudoJet> dijet;
        typedef std::pair<dijet, dijet>         hpair;

        vector<hpair> hcand = {hpair(dijet(srj[0], srj[1]), dijet(srj[2], srj[3])),
                               hpair(dijet(srj[0], srj[2]), dijet(srj[1], srj[3])),
                               hpair(dijet(srj[0], srj[3]), dijet(srj[1], srj[2]))};
        auto bc = std::min_element(hcand.begin(), hcand.end(), [](hpair& p1, hpair& p2) {
            return (fabs((p1.first.first + p1.first.second).m()
                         - (p1.second.first + p1.second.second).m())
                    < fabs((p2.first.first + p2.first.second).m()
                           - (p2.second.first + p2.second.second).m()));
        });

        // Construct pseudojets
        dijet&                   hp1     = (*bc).first;
        dijet&                   hp2     = (*bc).second;
        PseudoJet                higgs1  = hp1.first + hp1.second;
        PseudoJet                higgs2  = hp2.first + hp2.second;
        const fastjet::PseudoJet dihiggs = higgs1 + higgs2;

        // Sort by pT
        if (higgs1.pt() < higgs2.pt()) std::swap(higgs1, higgs2);
        if (hp1.first.pt() < hp1.second.pt()) std::swap(hp1.first, hp1.second);
        if (hp2.first.pt() < hp2.second.pt()) std::swap(hp2.first, hp2.second);

        // b-Tagging rate
        const int    nB = std::count(btags.begin(), btags.begin() + 4, BTAG); // Number of b's
        const int    nC = std::count(btags.begin(), btags.begin() + 4, CTAG); // Number of c's
        const int    nL = std::count(btags.begin(), btags.begin() + 4, LTAG); // Number of u,d,s,g
        const double selEff = btagProb(4, nB, nC, nL);
        const double selWgt = selEff * event_weight;
        // Reco fill
        HiggsFill(higgs1, higgs2, "res", 1, selWgt);

        const double diffHiggs_0 = fabs(higgs1.m() - 125.);
        const double diffHiggs_1 = fabs(higgs2.m() - 125.);
        const double massWindow  = 20.0;

        const bool signal0 = diffHiggs_0 < massWindow;
        const bool signal1 = diffHiggs_1 < massWindow;

        const bool control0 = !signal0 && diffHiggs_0 < 2.0 * massWindow;
        const bool control1 = !signal1 && diffHiggs_1 < 2.0 * massWindow;

        if (signal0 && signal1) {
            HiggsFill(higgs1, higgs2, "res", 2, selWgt);
            resNTuple << signal << "\t" << GetSample() << "\t" << selWgt << "\t" << higgs1.pt()
                      << "\t" << higgs2.pt() << "\t" << dihiggs.pt() << "\t" << higgs1.m() << "\t"
                      << higgs2.m() << "\t" << dihiggs.m() << "\t" << higgs1.delta_R(higgs2) << "\t"
                      << getDPhi(higgs1.phi(), higgs2.phi()) << "\t"
                      << fabs(higgs1.eta() - higgs2.eta()) << "\t" << Chi(higgs1, higgs2) << "\t"
                      << hp1.first.pt() << "\t" << hp1.second.pt() << "\t" << hp2.first.pt() << "\t"
                      << hp2.second.pt() << "\t" << std::endl;
        }
        else if ((signal0 && control1) || (signal1 && control0)) {
            HiggsFill(higgs1, higgs2, "res", 3, selWgt);
        }
        else if (control0 && control1) {
            HiggsFill(higgs1, higgs2, "res", 4, selWgt);
        }
        return selWgt;
    }

    return 0;
}

double OxfordSidebandAnalysis::IntermediateAnalysis(
      vector<PseudoJet> const& largeRJets, vector<PseudoJet> const& smallRJets,
      vector<vector<PseudoJet>> const& largeRsubJets, vector<vector<btagType>> const& largeRbtags,
      vector<btagType> const& smallRbtags, bool const& signal, double const& event_weight) {
    if (smallRJets.size() >= 2 && largeRJets.size() == 1 && largeRbtags[0].size() >= 2
        && largeRsubJets[0].size() >= 2) // MDT + reco cut
    {
        // Identify small-R jets separated from merged Higgs
        fastjet::Selector radialSelector = !SelectorCircle(1.2);
        radialSelector.set_reference(largeRJets[0]);
        const std::vector<fastjet::PseudoJet> separated =
              SelectorNHardest(2)(radialSelector(smallRJets));
        if (separated.size() != 2) return 0;

        // Construct the Higgs candidates
        fastjet::PseudoJet higgs1             = largeRJets[0];
        fastjet::PseudoJet higgs2             = separated[0] + separated[1];
        fastjet::PseudoJet res_lead_subjet    = separated[0];
        fastjet::PseudoJet res_sublead_subjet = separated[1];
        if (higgs1.pt() < higgs2.pt()) std::swap(higgs1, higgs2);
        if (res_lead_subjet.pt() < res_sublead_subjet.pt())
            std::swap(res_lead_subjet, res_sublead_subjet);

        // b-tagging weights
        const std::vector<btagType> separated_bTag = BTagging(separated);
        const int nB = std::count(largeRbtags[0].begin(), largeRbtags[0].begin() + 2, BTAG)
                       + (separated_bTag[0] == BTAG)
                       + (separated_bTag[1] == BTAG); // Number of true b-subjets
        const int nC = std::count(largeRbtags[0].begin(), largeRbtags[0].begin() + 2, CTAG)
                       + (separated_bTag[0] == CTAG)
                       + (separated_bTag[1] == CTAG); // Number of fake b-subjets
        const int nL = std::count(largeRbtags[0].begin(), largeRbtags[0].begin() + 2, LTAG)
                       + (separated_bTag[0] == LTAG)
                       + (separated_bTag[1] == LTAG); // Number of fake b-subjets
        const double selEff = btagProb(4, nB, nC, nL);
        const double selWgt = selEff * event_weight;

        if (nB + nC + nL != 4) return 0;

        HiggsFill(higgs1, higgs2, "inter", 1, selWgt);

        const double diffHiggs_0 = fabs(higgs1.m() - 125.);
        const double diffHiggs_1 = fabs(higgs2.m() - 125.);
        const double massWindow  = 20.0;

        const bool signal0 = diffHiggs_0 < massWindow;
        const bool signal1 = diffHiggs_1 < massWindow;

        const bool control0 = !signal0 && diffHiggs_0 < 2.0 * massWindow;
        const bool control1 = !signal1 && diffHiggs_1 < 2.0 * massWindow;

        if (signal0 && signal1) {
            HiggsFill(higgs1, higgs2, "inter", 2, selWgt);
            const PseudoJet dihiggs_inter = higgs1 + higgs2;
            const double    split12       = SplittingScales(largeRJets[0]);
            const double    ntau21        = tau21(largeRJets[0]);
            const double    nC2           = C2(largeRJets[0]);
            const double    nD2           = D2(largeRJets[0]);

            // Fill tuple
            intNTuple << signal << "\t" << GetSample() << "\t" << selWgt << "\t" << higgs1.pt()
                      << "\t" << higgs2.pt() << "\t" << dihiggs_inter.pt() << "\t" << higgs1.m()
                      << "\t" << higgs2.m() << "\t" << dihiggs_inter.m() << "\t"
                      << higgs1.delta_R(higgs2) << "\t" << getDPhi(higgs1.phi(), higgs2.phi())
                      << "\t" << fabs(higgs1.eta() - higgs2.eta()) << "\t" << Chi(higgs1, higgs2)
                      << "\t" << largeRsubJets[0][0].pt() << "\t" << largeRsubJets[0][1].pt()
                      << "\t" << res_lead_subjet.pt() << "\t" << res_sublead_subjet.pt() << "\t"
                      << split12 << "\t" << ntau21 << "\t" << nC2 << "\t" << nD2 << "\t"
                      << std::endl;
        }
        else if ((signal0 && control1) || (signal1 && control0)) {
            HiggsFill(higgs1, higgs2, "inter", 3, selWgt);
        }
        else if (control0 && control1) {
            HiggsFill(higgs1, higgs2, "inter", 4, selWgt);
        }

        return selWgt;
    }

    return 0;
}

// General fill for reconstructed higgs quantities
void OxfordSidebandAnalysis::HiggsFill(fastjet::PseudoJet const& H0, fastjet::PseudoJet const& H1,
                                       std::string const& analysis, size_t const& cut,
                                       double const& weight) {
    if (H0.pt() < H1.pt())
        std::cerr << "HiggsFill WARNING: pT ordering incorrect! " << analysis << "  " << cut << "  "
                  << H0.pt() << "  " << H1.pt() << std::endl;

    if (debug) std::cout << "HiggsFill INFO: cut = " << cut << std::endl;
    if (debug) std::cout << "HiggsFill INFO: analysis = " << analysis << std::endl;

    // Histo fill suffix
    const std::string suffix = "_" + analysis + cString[cut];

    // Record cutflow
    FillHistogram("CF_" + analysis, weight, cut + 0.1);
    FillHistogram("CFN_" + analysis, 1., cut + 0.1);

    // Histograms for reconstructed Higgs candidates
    FillHistogram("pt_H0" + suffix, weight, H0.pt());
    FillHistogram("pt_H1" + suffix, weight, H1.pt());

    FillHistogram("m_H0" + suffix, weight, H0.m());
    FillHistogram("m_H1" + suffix, weight, H1.m());

    FillHistogram("eta_H0" + suffix, weight, H0.eta());
    FillHistogram("eta_H1" + suffix, weight, H1.eta());

    FillHistogram("phi_H0" + suffix, weight, H0.phi());
    FillHistogram("phi_H1" + suffix, weight, H1.phi());

    FillHistogram("ptHptH" + suffix, weight, H0.pt(), H1.pt());
    FillHistogram("mHmH" + suffix, weight, H0.m(), H1.m());

    FillHistogram("dR_HH" + suffix, weight, H0.delta_R(H1));
    FillHistogram("dPhi_HH" + suffix, weight, getDPhi(H0.phi(), H1.phi()));
    FillHistogram("dEta_HH" + suffix, weight, fabs(H0.eta() - H1.eta()));
    FillHistogram("chi_HH" + suffix, weight, Chi(H0, H1));

    // Reconstruct di-Higgs system
    const fastjet::PseudoJet dihiggs = H0 + H1;

    // Histograms for reconstructed di-Higgs system
    FillHistogram("m_HH" + suffix, weight, dihiggs.m());
    FillHistogram("pt_HH" + suffix, weight, dihiggs.pt());
}

void OxfordSidebandAnalysis::BoostFill(fastjet::PseudoJet const& H, std::string const& analysis,
                                       size_t const& cut, double const& weight) {
    // Histo fill suffix
    const std::string suffix = "_" + analysis + cString[cut];

    // Splitting scales
    const double split12 = SplittingScales(H);

    // 2-subjettiness / 1-subjettiness
    const double btau21 = tau21(H);

    // C2/D2 energy correlation double-ratio
    const double bC2 = C2(H);
    const double bD2 = D2(H);

    FillHistogram("split12_fj" + suffix, weight, split12);
    FillHistogram("tau21_fj" + suffix, weight, btau21);
    FillHistogram("C2_fj" + suffix, weight, bC2);
    FillHistogram("D2_fj" + suffix, weight, bD2);
}

void OxfordSidebandAnalysis::BoostFill(fastjet::PseudoJet const& H0, fastjet::PseudoJet const& H1,
                                       std::string const& analysis, size_t const& cut,
                                       double const& weight) {
    if (H0.pt() < H1.pt())
        std::cerr << "BoostFill WARNING: pT ordering incorrect! " << analysis << "  " << cut << "  "
                  << H0.pt() << "  " << H1.pt() << std::endl;

    if (debug) std::cout << "BoostFill INFO: cut = " << cut << std::endl;
    if (debug) std::cout << "BoostFill INFO: analysis = " << analysis << std::endl;

    // Histo fill suffix
    const std::string suffix = "_" + analysis + cString[cut];

    // Splitting scales
    if (debug) std::cout << "BoostFill INFO: Calculate splitting scales" << std::endl;
    if (debug)
        std::cout << "BoostFill INFO! H0 pt = " << H0.pt() << " H1 pt = " << H1.pt() << std::endl;
    const double split12_fj1 = SplittingScales(H0);
    const double split12_fj2 = SplittingScales(H1);

    // 2-subjettiness / 1-subjettiness
    if (debug) std::cout << "BoostFill INFO: Calculate n-subjettiness" << std::endl;
    const double tau21_fj1 = tau21(H0);
    const double tau21_fj2 = tau21(H1);

    // C2 energy correlation double-ratio
    if (debug) std::cout << "BoostFill INFO: Calculate C2" << std::endl;
    const double C2_fj1 = C2(H0);
    const double C2_fj2 = C2(H1);

    // D2 energy correlation double-ratio
    if (debug) std::cout << "BoostFill INFO: Calculate D2" << std::endl;
    const double D2_fj1 = D2(H0);
    const double D2_fj2 = D2(H1);

    if (debug)
        std::cout << "BoostFill INFO! H0 split12 = " << split12_fj1
                  << " H1 split12 = " << split12_fj2 << std::endl;
    if (debug)
        std::cout << "BoostFill INFO! H0 tau21 = " << tau21_fj1 << " H1 tau21 = " << tau21_fj2
                  << std::endl;
    if (debug)
        std::cout << "BoostFill INFO! H0 C2 = " << C2_fj1 << " H1 C2 = " << C2_fj2 << std::endl;
    if (debug)
        std::cout << "BoostFill INFO! H0 D2 = " << D2_fj1 << " H1 D2 = " << D2_fj2 << std::endl;

    FillHistogram("split12_fj1" + suffix, weight, split12_fj1);
    FillHistogram("split12_fj2" + suffix, weight, split12_fj2);

    FillHistogram("tau21_fj1" + suffix, weight, tau21_fj1);
    FillHistogram("tau21_fj2" + suffix, weight, tau21_fj2);

    FillHistogram("C2_fj1" + suffix, weight, C2_fj1);
    FillHistogram("C2_fj2" + suffix, weight, C2_fj2);

    FillHistogram("D2_fj1" + suffix, weight, D2_fj1);
    FillHistogram("D2_fj2" + suffix, weight, D2_fj2);

    if (debug) std::cout << "BoostFill INFO: Finished BoostFill" << std::endl;
}
