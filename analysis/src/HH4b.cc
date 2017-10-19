/////////////////////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>
#include <functional>
#include <random>
#include <limits>
#include <locale>

#include "analysis.h"
#include "detector.h"
#include "run.h"
#include "samples.h"

using namespace Pythia8;

////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::uint32_t> initSeeds(runCard const& run, sampleCard const& sample,
                                     int const& subsample) {
    const uint32_t primary     = run.runseed;
    const uint32_t samplename  = std::hash<std::string>()(sample.samplename);
    const uint32_t runname     = std::hash<std::string>()(run.runname);
    const uint32_t subsampleID = subsample;

    std::seed_seq              seq{primary, samplename, runname, subsampleID};
    std::vector<std::uint32_t> seeds(4);
    seq.generate(seeds.begin(), seeds.end());
    return seeds;
}

////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
    // Set locale, so we have thousands separators
	std::locale my_locale{""};
	std::cout.imbue(my_locale);

    if (argc != 4) {
        cerr << "Error: Wrong number of arguments!" << endl;
        cerr << "Usage: HH4b <run card> <sample card> <subsample>" << endl;
        exit(-1);
    }

    // Read run card
    const std::string runfile = std::string(argv[1]);
    const runCard     run(runfile);

    // Read sample card
    const std::string samplefile = std::string(argv[2]);
    const sampleCard  sample(samplefile);

    // Determine subsample constants
    const int subsample   = atoi(argv[3]);
    const int sampleStart = subsample * run.sub_samplesize; // start point of the subsample

    const std::vector<std::uint32_t> seeds = initSeeds(run, sample, subsample);
    const uint32_t pythiaSeed = ((double)seeds[0] / pow(2, 32)) * 9E8; // Pythia seeds must be < 9E8
    const uint32_t pileupSeed = ((double)seeds[1] / pow(2, 32)) * 9E8; // Pythia seeds must be < 9E8
    const uint32_t detectorSeed = seeds[2];
    const uint32_t analysisSeed = seeds[3];

    int index_subsample_size = 0;
    std::vector<long long> subsample_indices;

    cout << "Processing sample: " << sample.samplename << ", subsample: " << subsample << std::endl;
    cout << "  RNG Seeds - Shower:   " << pythiaSeed << std::endl
         << "            - PU:       " << pileupSeed << std::endl
         << "            - Detector: " << detectorSeed << std::endl
         << "            - Analysis: " << analysisSeed << std::endl;

    // Initialise Pythia and HepMC
    Pythia        pythiaRun(std::string(PYTHIADIR));  // Pythia input
    std::ifstream hepmc_is(sample.eventpath.c_str()); // HepMC input

    // Initialise the event sample and weight normalisation
    double weight_norm = 0;
    if (!sample.hepmc)
        InitPythia(run, sample, pythiaSeed, pythiaRun, weight_norm);
    else
        InitHepMC(run, sample, weight_norm, index_subsample_size, subsample_indices);
    weight_norm *= 1000 * sample.xsec_norm; // Includes pb->fb conversion

    // Initialse Analyses and detector simulation
    vector<Analysis*> analyses;
    InitAnalyses(analyses, run, sample, subsample);
    Detector detector(run, sample, pileupSeed, detectorSeed);

    // Skip to subsample x
    if (!sample.hepmc) {
        cout << "Subsample " << subsample << ": skipping to event " << sampleStart << endl;
        for (int iEvent = 0; iEvent < sampleStart; ++iEvent) {
            double     dum;
            finalState dum2;
            get_final_state_particles(pythiaRun, dum2, dum);
        }
    }
    else {
	long long sampleStartByte = subsample_indices[subsample * (run.sub_samplesize / index_subsample_size)];
        cout << "Subsample " << subsample << ": fast-forwarding to event " << sampleStart
             << " (byte " << sampleStartByte << ")" << endl;
	if (subsample > 0) { // don't need to fast-forward if we're starting at the beginning
	    // HepMC *needs* to read the first event, so we do that, then rewind and skip lines
	    HepMC::GenEvent event;
	    cout << "Reading first event\n";
	    event.read(hepmc_is);
	    cout << "Fast-forwarding\n";
            hepmc_is.seekg(sampleStartByte);
	    cout << "Start analysis.\n";
	}
    }

    std::cout << "Number of analyses loaded: " << analyses.size() << std::endl;

    // Begin loop over events
    cout << "*************** Analysis Begins ***************" << endl;
    const int targetSize = min(run.sub_samplesize, sample.nevt_sample - sampleStart);
    cout << "Analysing: " << targetSize << " events" << endl;
    for (int iEvent = 0; iEvent < targetSize; ++iEvent) {
	bool extOutput = false;
        finalState ifs, fs; // The event final state

        double event_weight = 0;
        if (!sample.hepmc) {
            get_final_state_particles(pythiaRun, ifs, event_weight);
	}
        else {
            get_final_state_particles(hepmc_is, ifs, event_weight, extOutput);
	}

	if(extOutput) continue;
        // Perform detector simulation
        detector.Simulate(ifs, fs);
        // Run over analyses
        double gen_weight = event_weight; // To plot raw weights
        event_weight *= weight_norm;
        for (auto& analyse : analyses){
            analyse->Analyse(sample.is_signal, event_weight, fs, gen_weight);
        }

        if ((iEvent + 1) % 100 == 0) cout << iEvent + 1 << " events analysed" << endl;
    }

    // Clean up
    for (auto& analyse : analyses) delete analyse;
    hepmc_is.close();

    // End of the main program
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
