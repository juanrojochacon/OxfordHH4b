// newanalyis.h jr/nh 21/10/14
#pragma once

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <vector>

using std::string;

// Fwd fastjet
namespace fastjet {
class PseudoJet;
}

// Fwd YODA
namespace YODA {
class Histo1D;
class Histo2D;
}

class runCard;
class sampleCard;

// Final state particles
typedef std::vector<fastjet::PseudoJet> finalState;

class Analysis {
  public:
    Analysis(string const& name, runCard const&, sampleCard const&, int const& subsample);
    virtual ~Analysis();

    string const& GetName() const { return analysisName; };
    string const& GetRoot() const { return analysisRoot; };
    string const& GetSample() const { return sampleName; };
    int const&    GetSubSample() const { return subSample; }

    int const&    GetNPassed() const { return nPassed; };
    double const& GetWeight() const { return passedWeight; };

    double GetCutWeight() const { return cutWeight; };

    virtual void Analyse(bool const& signal, double const& weightnorm, finalState const&,
                         double gen_weight = 1.0) {
        totalWeight += weightnorm;
    };

    void Export();

    static bool Verbose;

  protected:
    void BookHistogram(YODA::Histo1D*, string const& name);
    void BookHistogram(YODA::Histo2D*, string const& name);

    void FillHistogram(string const& rname, double const& weight, double const& coord);
    void FillHistogram(string const& rname, double const& weight, double const& coord,
                       double const& coord2);

    void Cut(double const& weight) { cutWeight += weight; };
    void Pass(double const& weight) {
        nPassed++;
        passedWeight += weight;
    };

  private:
    // Standard analysis info
    const string analysisName; //!< Name of analysis
    const string analysisRoot; //!< Path to root of analysis results folder
    const string sampleName;   //!< Name of the current sample
    const int    subSample;    //!< Subsample name

    int    nPassed;      //!< Number of events passing analysis cuts
    double totalWeight;  //!< Total sample weight
    double passedWeight; //!< Total weight of passed events
    double cutWeight;    //!< Total weight explicitly cut

    std::map<int, YODA::Histo1D*> bookedHistograms_1D;
    std::map<int, YODA::Histo2D*> bookedHistograms_2D;
};
